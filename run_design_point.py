# -*- coding: utf-8 -*-
"""
설계점 (Design Point) 계산 스크립트 (run_design_point.py)

(수정됨: 2025-11-11 - Burner 아키텍처 수정)
(Tt4_R을 입력으로 사용하고, Burner.calculate_from_T4를 호출)
"""

import engine_components as eng
import cantera as ct

# --- 단위 변환 상수 ---
LBM_TO_KG = 0.453592
BTU_PER_LBM_TO_J_PER_KG = 2326.0

def print_state(name, state, W_pps=None):
    """엔진 스테이션의 상태를 Imperial 단위로 출력합니다."""
    T_R = state.T_K * eng.KELVIN_TO_RANKINE
    P_psia = state.P_Pa / eng.PSI_TO_PA
    h_btu_lbm = state.h_jkg * eng.J_PER_KG_TO_BTU_PER_LBM
    s_jkgK = state.s_jkgK # (SI 단위 유지)
    W_val = W_pps if W_pps is not None else state.W_pps
    
    print(f"--- {name} ---")
    print(f"  T: {T_R:,.2f} °R")
    print(f"  P: {P_psia:,.2f} psia")
    print(f"  h: {h_btu_lbm:,.2f} BTU/lbm")
    print(f"  s: {s_jkgK:,.2f} J/kg-K")
    print(f"  W: {W_val:,.2f} pps")


def run_design_point_check():
    """
    논문 Table 1 (이륙) 및 Table A1 (주기 데이터) 기반 검증
    """
    print("=== 설계점(Design Point) 계산 시작 ===")
    
    # 1. 입력값 정의 (Imperial → 즉시 SI 변환, 내부 계산은 SI로 일관 처리)
    # 환경 (Imperial)
    Alt_ft = 0.0
    MN_flight = 0.0
    T_amb_R = 518.67 # T_std
    P_amb_psia = 14.696 # P_std

    # 환경 (SI)
    T_amb_K = T_amb_R * eng.RANKINE_TO_KELVIN
    P_amb_Pa = P_amb_psia * eng.PSI_TO_PA
    
    # 컴포넌트 성능 (무차원/효율은 동일)
    PR_comp = 7.0
    Eff_comp = 0.87
    # 공기 질량 유량 (Imperial → SI)
    W_air_pps = 44.0
    W_air_kgps = W_air_pps * LBM_TO_KG
    
    # 연소기 (수정: Tt4_R을 입력으로 사용)
    Tt4_R_target = 2100.0  # 논문 설계점 목표 온도
    Burner_dP = 0.05 # 압력 손실 (비율)
    # LHV는 사용하지 않음 (참고용)
    LHV_btu_lbm = 18400.0
    
    # 터빈
    Eff_turb = 0.85
    
    # 노즐
    Cfg = 0.98 # 추력 계수
    Cd = 0.99 # 유량 계수 (추정치. 1.0에 가까움)
    
    # 2. Cantera 객체 설정
    air, fuel = eng.setup_gas()

    # 3. 계산 시작
    
    # --- 스테이션 0 (Ambient) ---
    air.TP = T_amb_K, P_amb_Pa
    
    state0 = eng.GasState(air, W_air_kgps)
    print_state("Station 0: Ambient", state0, W_pps=W_air_pps)

    # --- 스테이션 2 (Inlet) ---
    # 주의: Inlet.calculate_ram은 Imperial 입력을 기대하므로 경계에서만 변환
    Pt2_psia, Tt2_R = eng.Inlet.calculate_ram(P_amb_psia, T_amb_R, MN_flight, air)
    
    air.TP = Tt2_R * eng.RANKINE_TO_KELVIN, Pt2_psia * eng.PSI_TO_PA
    state2 = eng.GasState(air, W_air_kgps)
    print_state("Station 2: Inlet Face", state2, W_pps=W_air_pps)

    # --- 스테이션 3 (Compressor) ---
    Pt3_psia, Tt3_R, W_comp_jkg = eng.Compressor.calculate(
        Tt2_R, Pt2_psia, PR_comp, Eff_comp, air
    )
    
    air.TP = Tt3_R * eng.RANKINE_TO_KELVIN, Pt3_psia * eng.PSI_TO_PA
    state3 = eng.GasState(air, W_air_kgps)
    print_state("Station 3: Compressor Exit", state3, W_pps=W_air_pps)
    print(f"  Compressor Work: {W_comp_jkg * eng.J_PER_KG_TO_BTU_PER_LBM:,.2f} BTU/lbm")

    # --- 스테이션 4 (Burner) ---
    # (수정: calculate_from_T4 호출)
    # Tt4_R_target을 맞추기 위한 연료량(far)을 계산
    state4, far = eng.Burner.calculate_from_T4(
        state3, Tt4_R_target, Burner_dP, fuel
    )
    
    # 질량유량은 SI로 계산
    # W_fuel = W_air * far
    W_fuel_kgps = W_air_kgps * far
    W_total_kgps = W_air_kgps + W_fuel_kgps
    state4.W_kgps = W_total_kgps
    # 표시용(Imperial)
    W_fuel_pps = W_fuel_kgps / LBM_TO_KG
    
    Tt4_actual_R = state4.T_K * eng.KELVIN_TO_RANKINE
    print_state("Station 4: Turbine Inlet", state4)
    print(f"  Fuel/Air Ratio (far): {far:.4f}")
    print(f"  Fuel Flow (Wf): {W_fuel_pps:.4f} pps")
    print(f"  Actual Tt4: {Tt4_actual_R:,.2f} °R (Target: {Tt4_R_target:.1f} °R)")

    # --- 스테이션 5 (Turbine) ---
    # 압축기 소모일(W_comp_jkg)을 입력하여 일 평형을 맞춤
    state5 = eng.Turbine.calculate_from_work(
        state4, W_comp_jkg, Eff_turb, far
    )
    state5.W_kgps = state4.W_kgps
    print_state("Station 5: Turbine Exit", state5, W_pps=W_total_kgps / LBM_TO_KG)

    # --- 스테이션 8/9 (Nozzle) ---
    Fnet, Fg, Fd, A_th_calc = eng.Nozzle.calculate_design_point(
        state5,
        W_air_pps,         # 램 항력 계산용
        P_amb_psia,
        MN_flight,
        Cfg,
        Cd,
        T_amb_R
    )

    # --- 최종 성능 계산 ---
    # SFC (표시용) = (Wf_pps * 3600) / Fnet_lbf
    SFC = (W_fuel_pps * 3600) / Fnet
    
    print("\n=== 최종 성능 결과 (Design Point) ===")
    print(f"                        |  논문 (Table 1)  |  Cantera 모델  |")
    print(f"------------------------|------------------|----------------|")
    print(f"순 추력 (Fnet) [lbf]    |     2,850.0      |   {Fnet:,.1f}    |")
    print(f"연료 소비율 (SFC)       |       0.990      |     {SFC:.3f}      |")
    print(f"공기 유량 (Wa) [pps]    |      44.0        |     {W_air_pps:.1f}      |")
    print(f"압력비 (PR_comp)        |       7.0        |     {PR_comp:.1f}      |")
    print(f"터빈 입구 온도 (T4) [R] |     2,100.0      |   {Tt4_actual_R:,.1f}   |")
    print(f"  - 총 추력 (Fg) [lbf]   |        N/A       |   {Fg:,.1f}    |")
    print(f"  - 램 항력 (Fd) [lbf]   |        N/A       |     {Fd:,.1f}      |")
    print(f"  - 연료 유량 (Wf) [pps] |    ~0.79 (계산)  |     {W_fuel_pps:.4f}     |")
    
    print(f"\n--- 탈설계점(Off-Design)에 필요한 핵심 파라미터 ---")
    print(f"  계산된 노즐 목 면적 (A_th_calc) [m^2]: {A_th_calc:.4f}")
    print(f"  (이 값을 'run_off_design_solver.py'의 'A_th_m2_fixed' 변수에 복사하세요.)")


if __name__ == "__main__":
    run_design_point_check()