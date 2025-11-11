# -*- coding: utf-8 -*-
"""
엔진 컴포넌트 라이브러리 (engine_components.py)

(수정: 2025-11-11)
(1) Burner: 'calculate_from_far' (솔버용)와 'calculate_from_T4' (설계점용)로 분리.
(2) Inlet: .SH 속성 오류를 .setState_SH() 메서드로 수정.
(3) Turbine: .SH 속성 오류를 .setState_SH() 메서드로 수정하여 충돌 해결.
(4) Burner: J/kg -> J/kmol 단위 변환 오류 수정. (7004R 오류 해결)
"""

import cantera as ct
import numpy as np
from scipy.optimize import fsolve, brentq

# --- 단위 변환 상수 ---
PSI_TO_PA = 6894.76
RANKINE_TO_KELVIN = 1.0 / 1.8
KELVIN_TO_RANKINE = 1.8
LBF_TO_N = 4.44822
LBM_TO_KG = 0.453592
BTU_PER_LBM_TO_J_PER_KG = 2326.0
FT_PER_S_TO_M_PER_S = 0.3048
G_IMP = 32.174  # ft/s^2
G_SI = 9.80665 # m/s^2
J_PER_KG_TO_BTU_PER_LBM = 1.0 / 2326.0
P_STD_PA = 101325.0
T_STD_K = 288.15
P_STD_PSIA = P_STD_PA / PSI_TO_PA
T_STD_R = T_STD_K * KELVIN_TO_RANKINE


def setup_gas():
    """Cantera 가스 객체(공기, 연료)를 설정합니다."""
    # 1. 공기 객체 생성 (air.yaml 사용)
    air = ct.Solution('air.yaml')
    # Argon은 사용하지 않고 N2로 근사 (O2/N2만 사용)
    air.X = {'O2': 0.21, 'N2': 0.79}
    
    # 2. 연료 객체 생성 (nDodecane_Reitz.yaml 사용)
    # nDodecane_Reitz.yaml에는 두 개의 phase가 있음: nDodecane_RK, nDodecane_IG
    # ideal-gas phase를 사용
    fuel = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_IG')
    # n-dodecane (C12H26)을 Jet-A 연료의 대체물로 사용합니다.
    # 종 이름은 'c12h26' (소문자)
    fuel.X = {'c12h26': 1.0}
    fuel.TP = 298.15, ct.one_atm # 연료 주입 시 표준 상태
    return air, fuel

class GasState:
    """엔진 각 스테이션의 가스 상태를 저장하는 헬퍼 클래스"""
    def __init__(self, gas, W_kgps=0.0):
        self.gas = gas
        self.W_kgps = W_kgps # 질량 유량 (kg/s)
    
    @property
    def T_K(self): return self.gas.T
    @property
    def P_Pa(self): return self.gas.P
    @property
    def h_jkg(self): return self.gas.h # Cantera .h는 J/kg (질량 기준)
    @property
    def s_jkgK(self): return self.gas.s # Cantera .s는 J/kg-K (질량 기준)

    @property
    def T_R(self): return self.T_K * KELVIN_TO_RANKINE
    @property
    def P_psia(self): return self.P_Pa / PSI_TO_PA
    @property
    def h_btu_lbm(self): return self.h_jkg * J_PER_KG_TO_BTU_PER_LBM
    
    @property
    def W_pps(self): return self.W_kgps / LBM_TO_KG

    def update(self, T_K, P_Pa):
        """T, P로 상태 업데이트"""
        self.gas.TP = T_K, P_Pa

    def update_hP(self, h_jkg, P_Pa):
        """h(J/kg), P(Pa)로 상태 업데이트"""
        self.gas.HP = h_jkg, P_Pa

    def update_sP(self, s_jkgK, P_Pa):
        """s(J/kg-K), P(Pa)로 상태 업데이트"""
        self.gas.SP = s_jkgK, P_Pa


class Inlet:
    """
    스테이션 0 (Ambient) -> 2 (Inlet Face)
    """
    @staticmethod
    def calculate_ram(P_amb_psia, T_amb_R, mn, air):
        """
        램 회복(Ram Recovery)을 계산합니다.
        입력: Imperial (psia, R, mn), Cantera air 객체
        출력: Imperial (Pt, Tt)
        """
        if mn == 0:
            return P_amb_psia, T_amb_R

        T_amb_K = T_amb_R * RANKINE_TO_KELVIN
        P_amb_Pa = P_amb_psia * PSI_TO_PA
        
        # 1. 스테이션 0 (정적 상태) 설정
        air.TP = T_amb_K, P_amb_Pa
        h_amb_jkg = air.h
        s_amb_jkgK = air.s
        
        # 2. 비행 속도 계산 (V = mn * a)
        a_amb_mps = air.sound_speed # 음속
        V_flight_mps = mn * a_amb_mps
        
        # 3. 스테이션 2 (총량 상태) 계산
        # 총 엔탈피 = 정적 엔탈피 + 운동 에너지
        h_total_jkg = h_amb_jkg + 0.5 * V_flight_mps**2
        
        # 램 회복은 등엔트로피 과정으로 가정 (s_total = s_static)
        # (수정: .SH 속성 대신 .setState_SH() 메서드 사용)
        air.setState_SH(s_amb_jkgK, h_total_jkg)
        
        Pt_Pa = air.P
        Tt_K = air.T
        
        return Pt_Pa / PSI_TO_PA, Tt_K * KELVIN_TO_RANKINE


class Compressor:
    """스테이션 2 (Inlet) -> 3 (Compressor Exit)"""
    @staticmethod
    def calculate(Tt_in_R, Pt_in_psia, PR, eff, air):
        """
        압축기 계산
        입력: Imperial, 출력: Imperial (Pt_out, Tt_out), SI (W_comp_jkg)
        """
        Tt_in_K = Tt_in_R * RANKINE_TO_KELVIN
        Pt_in_Pa = Pt_in_psia * PSI_TO_PA
        air.TP = Tt_in_K, Pt_in_Pa
        h_in_jkg = air.h
        s_in_jkgK = air.s

        Pt_out_Pa = Pt_in_Pa * PR
        
        # 1. 이상적(등엔트로피) 출구 계산
        air.SP = s_in_jkgK, Pt_out_Pa
        h_out_ideal_jkg = air.h
        
        # 2. 실제 출구 계산 (효율 적용)
        W_comp_jkg = (h_out_ideal_jkg - h_in_jkg) / eff # W_comp_jkg는 J/kg (공기 1kg당)
        h_out_jkg = h_in_jkg + W_comp_jkg
        air.HP = h_out_jkg, Pt_out_Pa
        
        Tt_out_K = air.T
        
        return Pt_out_Pa / PSI_TO_PA, Tt_out_K * KELVIN_TO_RANKINE, W_comp_jkg


class Burner:
    """
    (수정됨: 2025-11-11 - Burner 함수 분리)
    스테이션 3 (Compressor Exit) -> 4 (Turbine Inlet)
    """
    
    @staticmethod
    def _perform_combustion(state3, far_target, dP_ratio, fuel):
        """
        (내부 헬퍼) far_target을 기반으로 HP 평형 연소를 수행하고,
        결과 'gas' 객체와 'far'를 반환합니다.
        """
        Pt_out_Pa = state3.P_Pa * (1.0 - dP_ratio)
        gas = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_IG')
        
        h_air_in_jkg = state3.h_jkg
        fuel.TP = 298.15, ct.one_atm 
        h_fuel_in_jkg = fuel.h # J/kg

        try:
            h_in_jkg = (h_air_in_jkg + far_target * h_fuel_in_jkg) / (1.0 + far_target)
            
            air_X = state3.gas.mole_fraction_dict()
            fuel_X = fuel.mole_fraction_dict()
            
            moles_air = 1.0 / state3.gas.mean_molecular_weight
            moles_fuel = far_target / fuel.mean_molecular_weight
            moles_total = moles_air + moles_fuel
            
            X_dict = {
                'o2': (moles_air * air_X.get('O2', 0)) / moles_total,
                'n2': (moles_air * air_X.get('N2', 0)) / moles_total,
                'c12h26': (moles_fuel * fuel_X.get('c12h26', 0)) / moles_total
            }
            
            gas.X = X_dict
            
            # (수정됨: 2025-11-11 - 치명적 단위 오류 수정)
            MW_unburned = gas.mean_molecular_weight # J/kg -> J/kmol 변환용
            H_in_jkmol = h_in_jkg * MW_unburned # J/kmol

            gas.HP = H_in_jkmol, Pt_out_Pa
            gas.equilibrate('HP')
            
            return gas
            
        except Exception as e:
            print(f"경고: _perform_combustion (far={far_target}) 계산 실패: {e}")
            gas.TPX = state3.T_K, Pt_out_Pa, state3.gas.X
            return gas # 실패 시 공기 상태 반환

    @staticmethod
    def calculate_from_far(state3, far_target, dP_ratio, fuel):
        """
        [탈설계점 솔버용]
        연소기 계산. far(연료-공기비)를 입력받아 T4를 *계산*합니다.
        """
        gas = Burner._perform_combustion(state3, far_target, dP_ratio, fuel)
        state4 = GasState(gas)
        return state4

    @staticmethod
    def calculate_from_T4(state3, Tt_out_R, dP_ratio, fuel):
        """
        [설계점 검증용]
        연소기 계산. T4 (Tt_out)를 맞추기 위한 연료량(far)을 솔버로 계산.
        """
        Tt_out_K = Tt_out_R * RANKINE_TO_KELVIN
        
        def _objective_func(far_guess):
            """(실제 연소 온도 - 목표 온도)를 반환하는 목적 함수"""
            if far_guess < 0.001 or far_guess > 0.1:
                return 1e6 
            
            gas = Burner._perform_combustion(state3, far_guess, dP_ratio, fuel)
            return gas.T - Tt_out_K

        far_guess_initial = 0.017 # 설계점 근처
        
        try:
            # fsolve 대신 더 안정적인 brentq 사용 (단조 증가 함수로 가정)
            far_actual = brentq(_objective_func, 0.001, 0.05, xtol=1e-6, rtol=1e-6, maxiter=100)
            ier = 1
        except Exception as e:
            # brentq가 실패하면 (예: Tt_out_K가 범위를 벗어남)
            ier = 0
            mesg = str(e)

        if ier != 1: 
            print(f"경고: Burner 솔버가 Tt4={Tt_out_R}R을 맞추는 데 실패했습니다. {mesg}")
            gas = Burner._perform_combustion(state3, far_guess_initial, dP_ratio, fuel)
            far_actual = far_guess_initial
        else:
            # 성공한 far로 최종 gas 객체 생성
            gas = Burner._perform_combustion(state3, far_actual, dP_ratio, fuel)
        
        state4 = GasState(gas)
        return state4, far_actual


class Turbine:
    """
    스테이션 4 (Turbine Inlet) -> 5 (Turbine Exit)
    """
    @staticmethod
    def calculate(state4, PR_turb, eff, far):
        """
        [탈설계점 솔버용] 터빈 팽창을 계산합니다.
        PR_turb (압력비)를 *입력*받아 W_turb_jkg (생산일)을 *반환*합니다.
        """
        # (수정) 원본 상태 보존을 위해 exhaust_gas 객체 복사
        exhaust_gas = ct.Solution(state4.gas.source, phase_id=state4.gas.name)
        exhaust_gas.TPX = state4.T_K, state4.P_Pa, state4.gas.X

        h_in_jkg = state4.h_jkg
        s_in_jkgK = state4.s_jkgK
        P_in_Pa = state4.P_Pa

        Pt_out_Pa = P_in_Pa / PR_turb
        
        # 1. 이상적(등엔트로피) 출구 계산
        exhaust_gas.SP = s_in_jkgK, Pt_out_Pa
        h_out_ideal_jkg = exhaust_gas.h
        
        # 2. 실제 출구 계산 (효율 적용)
        h_drop_ideal = h_in_jkg - h_out_ideal_jkg
        h_drop_real = h_drop_ideal * eff
        h_out_jkg = h_in_jkg - h_drop_real
        
        exhaust_gas.HP = h_out_jkg, Pt_out_Pa
        state5 = GasState(exhaust_gas)
        
        # 터빈 생산일 (J/kg, 배기가스 1kg당)
        W_turb_jkg = h_drop_real
        
        return state5, W_turb_jkg
    
    @staticmethod
    def calculate_from_work(state4, W_comp_jkg, eff, far):
        """
        (수정됨: 2025-11-11 - SH/SP 오류 수정)
        [설계점 검증용] 터빈 팽창을 계산합니다.
        W_comp_jkg (압축기 소모일)을 *입력*받아 일 평형을 맞추고,
        state5 (출구 상태, Pt5 포함)를 *반환*합니다.
        """
        # (수정) 원본 상태 보존을 위해 exhaust_gas 객체 복사
        exhaust_gas = ct.Solution(state4.gas.source, phase_id=state4.gas.name)
        exhaust_gas.TPX = state4.T_K, state4.P_Pa, state4.gas.X

        h_in_jkg = state4.h_jkg
        s_in_jkgK = state4.s_jkgK
        
        # 터빈이 생산해야 할 일 = 압축기 소모일
        # W_turb (J/kg_exhaust) * W_exhaust = W_comp (J/kg_air) * W_air
        # h_drop_turb_jkg (J/kg_exhaust) * (W_air * (1+far)) = W_comp_jkg * W_air
        h_drop_turb_jkg = W_comp_jkg / (1.0 + far)
        
        # 2. 실제 출구 엔탈피 (J/kg_exhaust)
        h_out_jkg = h_in_jkg - h_drop_turb_jkg
        
        # 3. 이상적 출구 엔탈피 (효율로 역산)
        h_out_ideal_jkg = h_in_jkg - (h_drop_turb_jkg / eff)
        
        # (수정됨: 2025-11-11 - AttributeError: .SH 수정)
        
        # 4. 이상적 팽창 후의 압력(Pt_out_Pa)을 찾기 위해 fsolve 사용
        #    목표: exhaust_gas.SP = s_in_jkgK, P_guess 일 때
        #          exhaust_gas.h == h_out_ideal_jkg 를 만족하는 P_guess를 찾는다.
        
        # (수정) Cantera 객체를 복사하여 솔버 전용으로 사용
        solver_gas = ct.Solution(state4.gas.source, phase_id=state4.gas.name)
        solver_gas.TPX = state4.T_K, state4.P_Pa, state4.gas.X
        
        def _find_P_from_SH(P_guess_Pa):
            """S, H(이상적)를 만족하는 P를 찾는 솔버 목적 함수"""
            if P_guess_Pa <= 0: # 음수 압력 방지
                return 1e6
            try:
                solver_gas.SP = s_in_jkgK, P_guess_Pa
                h_guess_jkg = solver_gas.h # J/kg
                error = h_guess_jkg - h_out_ideal_jkg
                return error
            except Exception:
                return 1e6 # 유효하지 않은 압력

        # 초기 P 추측값 (PR=3 근처)
        P_initial_guess = state4.P_Pa / 3.0 
        
        P_solution, infodict, ier, mesg = fsolve(
            _find_P_from_SH, 
            P_initial_guess, 
            full_output=True,
            xtol=1e-6
        )

        if ier != 1:
            print(f"경고: Turbine.calculate_from_work의 Pt5 솔버 실패: {mesg}")
            Pt_out_Pa = P_initial_guess
        else:
            Pt_out_Pa = P_solution[0]
            
        # 5. 실제 출구 상태 (h_out_jkg, Pt_out_Pa)로 T, s를 찾음
        exhaust_gas.HP = h_out_jkg, Pt_out_Pa
        state5 = GasState(exhaust_gas)
        
        return state5

class Nozzle:
    """
    스테이션 5 (Turbine Exit) -> 8/9 (Nozzle Exit)
    """
    
    @staticmethod
    def _calculate_throat(state_in, P_amb_Pa):
        """(내부 헬퍼) 노즐 목(throat)의 상태를 계산합니다."""
        # (수정) 원본 상태 보존을 위해 exhaust_gas 객체 복사
        exhaust_gas = ct.Solution(state_in.gas.source, state_in.gas.name)
        exhaust_gas.TPX = state_in.T_K, state_in.P_Pa, state_in.gas.X
        
        h_in_jkg = state_in.h_jkg
        s_in_jkgK = state_in.s_jkgK
        P_in_Pa = state_in.P_Pa
        
        # 목(throat)이 초크(choked)되었는지 확인
        exhaust_gas.SP = s_in_jkgK, P_in_Pa 
        gamma = exhaust_gas.cp_mass / exhaust_gas.cv_mass
        P_crit_Pa = P_in_Pa * (2 / (gamma + 1))**(gamma / (gamma - 1))
        
        if P_crit_Pa > P_amb_Pa:
            P_th_Pa = P_crit_Pa # Choked
        else:
            P_th_Pa = P_amb_Pa # Unchoked
        
        # 목(throat)에서의 상태 계산
        try:
            exhaust_gas.SP = s_in_jkgK, P_th_Pa
        except Exception as e:
            print(f"경고: Nozzle._calculate_throat SP 계산 실패 (P_th={P_th_Pa}, P_in={P_in_Pa}). {e}")
            # 비상시 unchoked 상태로 강제
            P_th_Pa = P_amb_Pa
            if P_th_Pa <= 0: P_th_Pa = 1.0 # 압력이 0이 되는 것 방지
            exhaust_gas.SP = s_in_jkgK, P_th_Pa

        h_th_jkg = exhaust_gas.h
        rho_th_kg_m3 = exhaust_gas.density
        
        h_drop = h_in_jkg - h_th_jkg
        if h_drop < 0: h_drop = 0 # 물리적 오류 방지
        
        V_th_mps = np.sqrt(2 * h_drop)
        
        return rho_th_kg_m3, V_th_mps, h_th_jkg
    
    @staticmethod
    def _calculate_thrust(state5, W_air_pps, P_amb_psia, mn_flight, Cfg, T_amb_R):
        """(내부 헬퍼) 추력을 계산합니다."""
        # (수정) 원본 상태 보존
        exhaust_gas = ct.Solution(state5.gas.source, phase_id=state5.gas.name)
        exhaust_gas.TPX = state5.T_K, state5.P_Pa, state5.gas.X
        
        h_in_jkg = state5.h_jkg
        s_in_jkgK = state5.s_jkgK
        W_total_kgps = state5.W_kgps 
        P_amb_Pa = P_amb_psia * PSI_TO_PA

        # --- 1. 추력 계산 (논문 Eq 13) ---
        try:
            exhaust_gas.SP = s_in_jkgK, P_amb_Pa
        except Exception as e:
             print(f"경고: Nozzle._calculate_thrust SP 계산 실패 (P_amb={P_amb_Pa}). {e}")
             exhaust_gas.TPX = state5.T_K, P_amb_Pa, state5.gas.X # 강제
        
        h_out_ideal_jkg = exhaust_gas.h
        h_drop_ideal = h_in_jkg - h_out_ideal_jkg
        if h_drop_ideal < 0: h_drop_ideal = 0 
            
        V_ideal_mps = np.sqrt(2 * h_drop_ideal)
        Fg_N = (W_total_kgps * V_ideal_mps) * Cfg
        
        # 램 항력 (Fd)
        air, _ = setup_gas()
        air.TP = T_amb_R * RANKINE_TO_KELVIN, P_amb_Pa
        a_amb_mps = air.sound_speed
        V_flight_mps = mn_flight * a_amb_mps
        W_air_kgps = W_air_pps * LBM_TO_KG
        Fd_N = W_air_kgps * V_flight_mps
        Fnet_N = Fg_N - Fd_N
        
        return Fnet_N / LBF_TO_N, Fg_N / LBF_TO_N, Fd_N / LBF_TO_N

    @staticmethod
    def calculate(state5, W_air_pps, P_amb_psia, mn_flight, Cfg, Cd, T_amb_R, A_th_m2):
        """
        [탈설계점 솔버용] 노즐 계산.
        A_th_m2 (고정된 목 면적)을 *입력*받아 W_calc (계산 유량)을 *반환*합니다.
        """
        # 1. 추력 계산
        Fnet_lbf, Fg_lbf, Fd_lbf = Nozzle._calculate_thrust(
            state5, W_air_pps, P_amb_psia, mn_flight, Cfg, T_amb_R
        )
        
        # 2. 이론적 질량 유량 계산 (W_calc)
        P_amb_Pa = P_amb_psia * PSI_TO_PA
        rho_th_kg_m3, V_th_mps, _ = Nozzle._calculate_throat(state5, P_amb_Pa)
        
        # W_calc = (rho * V * A_th) * Cd
        W_calc_kgps = (rho_th_kg_m3 * V_th_mps * A_th_m2) * Cd
        
        return Fnet_lbf, Fg_lbf, Fd_lbf, W_calc_kgps / LBM_TO_KG

    @staticmethod
    def calculate_design_point(state5, W_air_pps, P_amb_psia, mn_flight, Cfg, Cd, T_amb_R):
        """
        [설계점 검증용] 노즐 계산.
        A_th_m2 (목 면적)을 *반환*합니다.
        """
        # 1. 추력 계산
        Fnet_lbf, Fg_lbf, Fd_lbf = Nozzle._calculate_thrust(
            state5, W_air_pps, P_amb_psia, mn_flight, Cfg, T_amb_R
        )
        
        # 2. 목 면적(A_th) 역산
        P_amb_Pa = P_amb_psia * PSI_TO_PA
        rho_th_kg_m3, V_th_mps, _ = Nozzle._calculate_throat(state5, P_amb_Pa)
        
        # W_total = A_th * (rho_th * V_th) * Cd
        W_total_kgps = state5.W_kgps
        if (rho_th_kg_m3 * V_th_mps) == 0:
             raise ValueError("Nozzle.calculate_design_point: Throat velocity or density is zero.")
             
        A_th_m2_calc = W_total_kgps / (rho_th_kg_m3 * V_th_mps * Cd)
        
        return Fnet_lbf, Fg_lbf, Fd_lbf, A_th_m2_calc