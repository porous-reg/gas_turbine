# -*- coding: utf-8 -*-
"""
탈설계점(Off-Design) 솔버 (run_off_design_solver.py)

(수정됨: 2-변수 비선형 솔버 구현)
(수정: 2025-11-11 - Burner.calculate_from_far 호출)
(수정: 2025-11-11 - 최종 성능 출력 로직 추가)
"""

import numpy as np
from scipy.optimize import root
import engine_components as eng
from engine_components import (
    Inlet, Compressor, Burner, Turbine, Nozzle,
    GasState, setup_gas,
    LBM_TO_KG, LBF_TO_N, RANKINE_TO_KELVIN, KELVIN_TO_RANKINE,
    PSI_TO_PA, P_STD_PSIA, T_STD_R, P_STD_PA, T_STD_K
)

# -------------------------------------------------------------
# TODO: (필수) 실제 성능 맵 데이터 구현
# -------------------------------------------------------------
# J85 압축기/터빈 맵 데이터가 필요합니다.
# 현재는 솔버 구조를 테스트하기 위한 임시 '자리 표시자'입니다.

def get_compressor_map(Nc_rpm, R_line):
    """
    (임시) 압축기 맵. Nc, R_line -> PR, Wc, Eff
    Nc_rpm: 실제 축 속도 (예: 16540)
    R_line: 작동 라인 (솔버가 찾는 미지수, 0~1.5)
    """
    # Nc_norm (0~1.02)
    Nc_norm = Nc_rpm / 16540.0 
    
    # 이 부분은 실제 2D 보간(interp2d)으로 대체되어야 합니다.
    PR_comp = (3.0 + 4.0 * Nc_norm) * R_line
    Wc_pps_norm = (15 + 29 * Nc_norm) * (2.0 - R_line) # R_line이 1일때 44
    Eff_comp = 0.87 - (1.0 - Nc_norm)**2 - (1.0 - R_line)**2 # 1,1에서 최대
    
    # 값 보장
    if PR_comp < 1.0: PR_comp = 1.0
    if Wc_pps_norm < 1.0: Wc_pps_norm = 1.0
    if Eff_comp < 0.5: Eff_comp = 0.5
    
    return PR_comp, Wc_pps_norm, Eff_comp

def get_turbine_map(Nc_rpm, R_line):
    """
    (임시) 터빈 맵. PR_turb, Eff_turb 반환
    (터빈 맵은 복잡하지만, 여기서는 압축기와 연동된다고 가정)
    """
    Nc_norm = Nc_rpm / 16540.0
    
    PR_turb = (1.5 + 1.163 * Nc_norm) * (0.8 + R_line * 0.2) # R_line 1에서 2.663
    Eff_turb = 0.85 - (1.0 - Nc_norm)**2
    
    return PR_turb, Eff_turb
# -------------------------------------------------------------

def uncorrect_flow(Wc_pps, Pt_in_psia, Tt_in_R):
    """ Wc (보정 유량) -> W (실제 유량) 변환 """
    theta = Tt_in_R / T_STD_R
    delta = Pt_in_psia / P_STD_PSIA
    W_pps = Wc_pps * delta / np.sqrt(theta)
    return W_pps

def calculate_performance(params, A_th_m2_fixed, x):
    """
    솔버가 호출할 목적 함수(Objective Function).
    미지수 x = [Nc_rpm, R_line] 를 입력받아,
    오차 [err_work, err_flow] 와 {performance_data}를 반환합니다.
    """
    Nc_rpm_guess, R_line_guess = x
    
    try:
        # --- 1. 가스 객체 설정 ---
        air, fuel = setup_gas()
        
        # --- 2. 고정 입력 (환경) ---
        P_amb_psia = params['P_amb_psia']
        T_amb_R = params['T_amb_R']
        mn_flight = params['mn_flight']
        
        # --- 3. Inlet (스테이션 0 -> 2) ---
        Pt2_psia, Tt2_R = Inlet.calculate_ram(P_amb_psia, T_amb_R, mn_flight, air)
        
        # --- 4. 압축기 (스테이션 2 -> 3) ---
        # 맵 조회 (미지수 사용)
        PR_comp, Wc_pps_norm, Eff_comp = get_compressor_map(Nc_rpm_guess, R_line_guess)
        
        # 실제 공기 유량 계산
        W_air_pps = uncorrect_flow(Wc_pps_norm, Pt2_psia, Tt2_R)
        W_air_kgps = W_air_pps * LBM_TO_KG
        
        # 압축기 계산
        Pt3_psia, Tt3_R, W_comp_jkg = Compressor.calculate(
            Tt2_R, Pt2_psia, PR_comp, Eff_comp, air
        )
        # 압축기 총 소모 일 (Power)
        Power_comp_W = W_comp_jkg * W_air_kgps
        
        # --- 5. 연소기 (스테이션 3 -> 4) ---
        state3 = GasState(air, W_air_kgps)
        state3.update(Tt3_R * RANKINE_TO_KELVIN, Pt3_psia * PSI_TO_PA)
        
        # (수정됨: calculate_from_far 호출)
        # 'far_target'을 직접 입력으로 사용
        state4 = Burner.calculate_from_far(
            state3, params['far_target'], params['Burner_dP'], fuel
        )
        far = params['far_target']
        
        # --- 6. 터빈 (스테이션 4 -> 5) ---
        # 맵 조회 (미지수 사용)
        PR_turb, Eff_turb = get_turbine_map(Nc_rpm_guess, R_line_guess)
        
        # 터빈 계산
        state5, W_turb_jkg = Turbine.calculate(
            state4, PR_turb, Eff_turb, far
        )
        
        # 터빈 총 생산 일 (Power)
        W_total_kgps = W_air_kgps * (1.0 + far)
        Power_turb_W = W_turb_jkg * W_total_kgps
        
        # --- 7. 노즐 (스테이션 5 -> 9) ---
        state5.W_kgps = W_total_kgps # 유량 설정
        
        Fnet_lbf, Fg_lbf, Fd_lbf, W_calc_pps = Nozzle.calculate(
            state5, W_air_pps, P_amb_psia, mn_flight, 
            params['Cfg'], params['Cd'], T_amb_R, 
            A_th_m2_fixed # 고정된 노즐 면적 사용
        )
        
        # --- 8. 오차 계산 ---
        err_work = (Power_turb_W - Power_comp_W) / 1e5 # 일 평형 오차 (스케일링)
        err_flow = (W_air_pps - W_calc_pps) / 100.0   # 유량 평형 오차 (스케일링)
        
        # 디버깅 출력 (필요시 주석 해제)
        # print(f"Guess x=({Nc_rpm_guess:.1f}, {R_line_guess:.3f}) -> Errors=({err_work:.4e}, {err_flow:.4e})")

        # (수정됨) 계산된 성능 데이터 반환
        W_fuel_pps = W_air_pps * far
        SFC = (W_fuel_pps * 3600) / Fnet_lbf if Fnet_lbf > 0 else 0.0

        performance_data = {
            "Fnet_lbf": Fnet_lbf,
            "SFC": SFC,
            "W_air_pps": W_air_pps,
            "W_fuel_pps": W_fuel_pps,
            "far_target": far,
            "Nc_rpm": Nc_rpm_guess,
            "R_line": R_line_guess,
            "PR_comp": PR_comp,
            "PR_turb": PR_turb,
            "Tt4": state4.T_R,
            "Pt4": state4.P_psia,
            "Pt5": state5.P_psia, "Tt5": state5.T_R,
        }

        return [err_work, err_flow], performance_data
        
    except Exception as e:
        # Cantera 오류 등 계산 실패 시 큰 오차 반환
        print(f"솔버 오류 발생: {e} (x={x})")
        return [1e6, 1e6], {}


def find_operating_point():
    """
    비선형 솔버를 설정하고 실행하여 작동점을 찾습니다.
    """
    
    # --- 고정 입력 파라미터 ---
    params = {
        'P_amb_psia': P_STD_PSIA, # 0 ft
        'T_amb_R': T_STD_R,       # 0 ft
        'mn_flight': 0.0,         # SLS
        'far_target': 0.0171,     # 스로틀 입력 (이륙 조건 far 근사치)
        'Burner_dP': 0.05,
        'Cfg': 0.98,
        'Cd': 1.0,
    }
    
    # 고정된 노즐 목 면적 (m^2)
    # (주의: 이 값은 run_design_point.py에서 정확히 역산된 값이어야 함)
    # J85-GE-13 근사치 사용 (0.088 m^2)
    A_th_m2_fixed = 0.0882 # (예시: run_design_point.py에서 0.0882가 나왔다고 가정)

    # --- 솔버 설정 ---
    # 미지수 x = [Nc_rpm, R_line]
    # 초기 추측값 (이륙 조건 근처)
    x0 = [16540.0, 1.0] 
    
    print(f"솔버 시작... 초기 추측값 x0 = {x0}")
    print(f"입력 조건: Alt=0ft, MN=0, FAR={params['far_target']:.5f}")
    
    # (수정됨) 솔버에 전달할 목적 함수를 래핑(wrapping)
    def objective_function(x):
        errors, _ = calculate_performance(params, A_th_m2_fixed, x)
        return errors
        
    # 솔버 실행
    # 'objective_function'이 [0, 0]이 되는 x를 찾습니다.
    sol = root(objective_function, x0, method='lm')
    
    # --- 결과 출력 ---
    if sol.success:
        print("\n--- 솔버 성공 ---")
        Nc_rpm_final, R_line_final = sol.x
        print(f"작동점 발견:")
        print(f"  - Nc (축 속도): {Nc_rpm_final:.2f} RPM")
        print(f"  - R_line (작동선): {R_line_final:.4f}")
        
        # (수정됨: 최종 성능 출력 로직 추가)
        # 최종 x값을 사용하여 전체 성능(Fnet 등)을 계산하고 출력
        errors, final_perf = calculate_performance(params, A_th_m2_fixed, sol.x)
        
        print("\n--- 최종 성능 (Off-Design) ---")
        print(f"  순 추력 (Fnet): {final_perf['Fnet_lbf']:,.2f} lbf")
        print(f"  연료 소비율 (SFC): {final_perf['SFC']:.4f}")
        print(f"  공기 유량 (Wa): {final_perf['W_air_pps']:.2f} pps")
        print(f"  연료 유량 (Wf): {final_perf['W_fuel_pps']:.4f} pps (far={final_perf['far_target']:.4f})")
        print(f"  압축기 PR: {final_perf['PR_comp']:.3f}")
        print(f"  터빈 PR: {final_perf['PR_turb']:.3f}")
        print("\n--- 스테이션 데이터 (Imperial) ---")
        print(f"  St 4 (Burner): P={final_perf['Pt4']:.2f} psia, T={final_perf['Tt4']:.2f} R")
        print(f"  St 5 (Turb):   P={final_perf['Pt5']:.2f} psia, T={final_perf['Tt5']:.2f} R")
        
    else:
        print("\n--- 솔버 실패 ---")
        print(f"메시지: {sol.message}")
        print("참고: 임시(Placeholder) 맵 데이터가 부정확하여 솔버가 실패할 수 있습니다.")
        print("      실제 맵 데이터를 'get_compressor_map'에 구현해야 합니다.")

if __name__ == "__main__":
    find_operating_point()