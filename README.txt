Python / Cantera 기반 고충실도 터보젯 모델

이 프로젝트는 NASA/TM-2017-219520 논문을 기반으로 J85 터보젯 엔진의 성능을 Python과 Cantera를 사용해 모델링합니다.

이 모델은 공학적 검증 및 머신러닝(ML) 데이터 생성을 목적으로 하는 고충실도(high-fidelity) 백엔드 프로그램입니다.

핵심 로직: 2가지 실행 모드

이 프로그램은 두 가지 서로 다른 모드로 작동합니다.

모드 1: 설계점 검증 (run_design_point.py)

목적: 우리 모델(engine_components.py)이 논문에 명시된 이륙 조건(Table 1) 데이터와 일치하는지 **검증(Validate)**합니다.

로직: 논문에 주어진 PR_comp, Eff_comp, W_air_pps 등 모든 값을 입력으로 사용합니다. Work_Compressor = Work_Turbine을 가정하여 터빈 출구 압력(Pt5)을 역산합니다.

핵심 결과: 논문과 일치하는 Fnet(추력), SFC(연료소비율)를 출력하고, 이 엔진의 고유한 물리적 특성인 **A_th_calc(노즐 목 면적)**을 계산해냅니다.

모드 2: 탈설계점 계산 (run_off_design_solver.py)

목적: 설계점 외의 모든 조건(다른 고도, 마하, 스로틀)에서의 엔진 성능을 **계산(Calculate)**합니다. (ML 데이터 생성용)

로직: scipy.optimize.root 비선형 솔버를 사용하여 물리 법칙을 만족하는 작동점(Operating Point)을 찾습니다.

핵심 입력:

고정 입력: Alt (고도), MN (마하 수), T4_R (터빈 입구 온도, 스로틀에 해당), 그리고 모드 1에서 계산한 A_th_m2_fixed (노즐 목 면적).

맵 데이터 (TODO): get_compressor_map, get_turbine_map 함수에 실제 맵 데이터 구현이 필요합니다.

솔버의 목표: 다음 두 가지 오차(error)가 0이 되는 미지수 x = [Nc_rpm, R_line] (축 속도, 압축기 작동선)를 찾습니다.

err_work = Power_Turbine - Power_Compressor (일 평형: 터빈이 만든 일과 압축기가 소모한 일이 같아야 함)

err_flow = W_in_compressor - W_calc_nozzle (유량 평형: 압축기로 들어간 공기 유량이 고정된 노즐 면적으로 나갈 수 있는 유량과 같아야 함)

핵심 결과: 솔버가 찾은 Nc_rpm (축 속도) 및 이 조건에서의 Fnet, SFC 등 최종 성능입니다.

파일 구성

engine_components.py:

로직: 엔진의 각 구성 요소(Inlet, Compressor, Burner, Turbine, Nozzle)를 Python 클래스/함수로 정의한 핵심 라이브러리입니다. Cantera를 사용하여 정밀한 열역학 계산을 수행합니다.

특징: Turbine과 Nozzle 클래스는 두 가지 모드를 모두 지원하기 위해 calculate() (솔버용)와 calculate_from_work() / calculate_design_point() (설계점 검증용) 함수를 별도로 가지고 있습니다.

run_design_point.py:

로직: 모드 1을 실행합니다.

필요한 입력: (코드 내 하드코딩) 논문 Table 1의 이륙 조건 값 (PR_comp = 7.0, Eff_comp = 0.87, W_air_pps = 44.0, Tt4_R = 2100.0 등)

결과 (출력): 콘솔에 논문 값과 비교된 Fnet, SFC를 출력하고, A_th_calc (노즐 목 면적) 값을 출력합니다.

run_off_design_solver.py:

로직: 모드 2를 실행합니다.

필요한 입력:

(코드 내 하드코딩) params 딕셔너리의 환경 변수 (Alt, MN, T4_R_target).

(사용자 수정 필요) A_th_m2_fixed: run_design_point.py 실행 결과로 나온 A_th_calc 값을 여기에 복사/붙여넣기 해야 합니다.

(TODO) get_compressor_map, get_turbine_map 함수 내에 실제 맵 데이터와 보간 로직을 구현해야 합니다. (현재는 임시 값 사용)

결과 (출력): 콘솔에 솔버가 찾은 Nc_rpm (축 속도)과 R_line (작동선)을 출력합니다.

실행 워크플로우 (Workflow)

1단계: 환경 설정 (최초 1회)

이 모델은 Python 3와 Cantera, Numpy, Scipy 라이브러리가 필요합니다.

# Cantera 설치 (Conda 사용 권장)
conda install -c cantera cantera

# Scipy, Numpy 설치 (Conda에 포함되어 있거나 pip로 설치)
pip install scipy numpy


2단계: 설계점 검증 및 노즐 면적(A_th) 확보

run_design_point.py를 실행하여 모델을 검증하고 A_th 값을 얻습니다.

python run_design_point.py


출력 예시:

=== 최종 성능 결과 (Design Point) ===
                        |  논문 (Table 1)  |  Cantera 모델  |
------------------------|------------------|----------------|
순 추력 (Fnet) [lbf]    |     2,850.0      |   2,851.5    |
연료 소비율 (SFC)       |       0.990      |     0.989      |
...
--- 탈설계점(Off-Design)에 필요한 핵심 파라미터 ---
  계산된 노즐 목 면적 (A_th_calc) [m^2]: 0.0882
  (이 값을 'run_off_design_solver.py'의 'A_th_m2_fixed' 변수에 복사하세요.)


3단계: 탈설계점(Off-Design) 솔버 설정

run_off_design_solver.py 파일을 열고, A_th_m2_fixed 변수 값을 2단계에서 얻은 값으로 수정합니다.

수정 전:
A_th_m2_fixed = 0.088

수정 후 (예시):
A_th_m2_fixed = 0.0882

4단계: (TODO) 성능 맵 데이터 구현

이 단계는 고충실도 계산을 위해 필수적입니다.

J85 압축기 및 터빈 맵 데이터(논문 Fig 10, 15 등)를 .csv 파일 등으로 구합니다.

run_off_design_solver.py의 get_compressor_map, get_turbine_map 함수 내부를 scipy.interpolate.interp2d (또는 RegularGridInterpolator)를 사용하여 실제 맵 데이터를 보간하도록 수정합니다.

5단계: 탈설계점(Off-Design) 솔버 실행

run_off_design_solver.py를 실행하여 params 딕셔너리에 설정된 조건(예: 30,000 ft, 0.8 MN, 90% 스로틀)의 성능을 계산합니다.

python run_off_design_solver.py


솔버가 성공하면(참고: 현재는 임시 맵 데이터로 인해 실패할 수 있음), 해당 조건의 Nc_rpm (축 속도)이 출력됩니다.