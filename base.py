# 전처리
import pandas as pd
train_df = pd.read_csv('train.csv')
## gdc 데이터 이용
import os
import pandas as pd
import glob

# 고유 SUBCLASS 이름 가져오기
cancer = train_df['SUBCLASS'].unique()

# 'SUBCLASS'와 'ID'를 제외한 열들 가져오기
remaining_columns = train_df.drop(columns=['SUBCLASS', 'ID']).columns

# 전체 데이터를 합치기 위한 리스트
all_dataframes = []

for item in cancer:
    # CSV 파일 경로 설정 (이전에 생성된 CSV 파일들)
    csv_file = f'{item}.csv'

    # 파일이 존재하는지 확인
    if not os.path.exists(csv_file):
        print(f"Skipping {csv_file}, as it does not exist.")
        continue

    # CSV 파일 불러오기
    df = pd.read_csv(csv_file)
    print(f"Loaded {csv_file} with shape {df.shape}")

    # 'Hugo_Symbol'이 remaining_columns와 일치하는 행만 남김
    df_filtered = df[df['Hugo_Symbol'].isin(remaining_columns)]

    # 필요한 열만 선택 (Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)
    df_filtered = df_filtered[['Hugo_Symbol', 'HGVSp_Short', 'Tumor_Sample_Barcode']]

    df_filtered['HGVSp_Short'] = df_filtered['HGVSp_Short'].str.strip()

    # 'HGVSp_Short'에서 'p.' 제거
    df_filtered['HGVSp_Short'] = df_filtered['HGVSp_Short'].str.replace(r'p\.', '', regex=False)

    # '='가 포함된 경우, 알파벳과 숫자를 분리하고 '='를 알파벳으로 바꾸기
    df_filtered['HGVSp_Short'] = df_filtered['HGVSp_Short'].str.replace(r'([A-Z])(\d+)=', r'\1\2\1', regex=True)

    # '='가 포함되거나 공백인 행 제거
    df_cleaned = df_filtered[~df_filtered['HGVSp_Short'].str.contains('=', na=False) & df_filtered['HGVSp_Short'].str.strip().ne('')]
    
    df_cleaned = df_cleaned[~df_cleaned['HGVSp_Short'].str.contains(r'\?', na=False)]
    # 새로운 열 'SUBCLASS' 추가하고 각 행에 item 값을 넣음
    df_cleaned['SUBCLASS'] = item
    
    # 처리된 데이터프레임을 리스트에 추가
    all_dataframes.append(df_cleaned)
    print(f"Cleaned {csv_file}, new shape: {df_cleaned.shape}")

# 모든 데이터를 하나로 합치기
if all_dataframes:
    ALL_DATA = pd.concat(all_dataframes, ignore_index=True)
    ALL_DATA.to_csv('ALL_DATA.csv', index=False)
    print(f"All files combined and saved to 'ALL_DATA.csv' with shape {ALL_DATA.shape}")
else:
    print("No dataframes to combine.")

# 남아있는 '='가 있는 행 확인
remaining_equals = ALL_DATA[ALL_DATA['HGVSp_Short'].str.contains('=', na=False)]
if not remaining_equals.empty:
    print("Rows with '=' remaining:")
    print(remaining_equals)
else:
    print("No rows with '=' remaining.")

a_empty_rows = ALL_DATA[ALL_DATA['Tumor_Sample_Barcode'].isna()]
a_empty_rows
# '?' 문자가 포함된 행 제거
ALL_DATA = ALL_DATA[~ALL_DATA['HGVSp_Short'].isna()]
ALL_DATA
#591756
# 'HGVSp_Short' 열에서 'fs' 뒤의 부분을 제거
ALL_DATA['HGVSp_Short'] = ALL_DATA['HGVSp_Short'].str.replace(r'fs.*', 'fs', regex=True)


# 변경된 결과 확인
print(ALL_DATA)

# 'HGVSp_Short' 열에서 'p.' 제거
ALL_DATA['HGVSp_Short'] = ALL_DATA['HGVSp_Short'].str.strip().str.replace('p.', '', regex=False)
# 그룹화하여 Protein_Change를 합치기
ALL_DATA = ALL_DATA.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode', 'SUBCLASS'], as_index=False).agg({
    'HGVSp_Short': lambda x: ' '.join(x.dropna().astype(str))  # NaN 값 제거 및 문자열로 변환 후 합치기
})

# 결과 출력
ALL_DATA
#502228
# 1. Hugo_Symbol 값의 고유 리스트 생성
hugo_symbols = train_df.drop(columns=['ID', 'SUBCLASS']).columns.tolist()

# 2. 새로운 DataFrame 초기화
result_rows = []  # 빈 리스트 초기화

# 3. Tumor_Sample_Barcode 및 SUBCLASS별로 데이터 처리
for (barcode, subclass), group in final_df.groupby(['Tumor_Sample_Barcode', 'SUBCLASS']):
    # 기본값을 가지는 행 생성 (모든 Hugo_Symbol은 'WT'로 설정)
    result_row = {hugo: 'WT' for hugo in hugo_symbols}
    
    # SUBCLASS와 ID 열 추가
    result_row['SUBCLASS'] = subclass
    result_row['ID'] = barcode  # Tumor_Sample_Barcode를 ID로 사용
    
    # 그룹 내의 각 Hugo_Symbol에 대해 Protein Change 값을 설정
    for _, row in group.iterrows():
        hugo_symbol = row['Hugo_Symbol']
        protein_change = row['HGVSp_Short']
        
        # 값이 이미 존재하는 경우, 띄어쓰기를 추가하여 값을 연결
        if result_row[hugo_symbol] != 'WT':
            result_row[hugo_symbol] += ' ' + protein_change  # 기존 값에 띄어쓰기를 추가하여 연결
        else:
            result_row[hugo_symbol] = protein_change  # 기본값을 설정
            
    result_rows.append(result_row)  # 결과 행을 리스트에 추가

# 4. 새로운 DataFrame 생성
result_df = pd.DataFrame(result_rows)

# 5. 빈 값은 'WT'로 대체
result_df.replace('', 'WT', inplace=True)

# 6. train_df의 열 추가 (없는 열을 'WT'로 채움)
for column in train_df.columns:
    if column not in result_df.columns:
        result_df[column] = 'WT'

result_df.to_csv("ALL_DATA_FINAL.csv",index=False)
# phantasus 데이터 이용
import pandas as pd
import requests

cancer = train_df['SUBCLASS'].unique()
# 데이터프레임을 저장할 리스트 초기화
dfs = []

# cancer 리스트에 있는 각 아이템에 대해 반복
for item in cancer:
    # 파일 경로를 지정
    file_path = f"mutations_merged_{item}.maf.txt"
    
    # 파일을 DataFrame으로 읽기 (탭으로 구분된 파일)
    df = pd.read_csv(file_path, delimiter='\t')
    
    # 필요한 열만 선택
    df = df[['Hugo_Symbol', 'Protein_Change','Tumor_Sample_Barcode']]
    df['SUBCLASS'] = item
    # 리스트에 DataFrame 추가
    dfs.append(df)

# 모든 DataFrame을 하나로 이어 붙이기
final_df = pd.concat(dfs, ignore_index=True)
a_empty_rows = final_df[final_df['Tumor_Sample_Barcode'].isna()]
a_empty_rows
# 필요한 열만 선택 (Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)
final_df = final_df[['Hugo_Symbol', 'Protein_Change', 'Tumor_Sample_Barcode','SUBCLASS']]

final_df['Protein_Change'] = final_df['Protein_Change'].str.strip()

# 'HGVSp_Short'에서 'p.' 제거
final_df['Protein_Change'] = final_df['Protein_Change'].str.replace(r'p\.', '', regex=False)

# '='가 포함된 경우, 알파벳과 숫자를 분리하고 '='를 알파벳으로 바꾸기
final_df['Protein_Change'] = final_df['Protein_Change'].str.replace(r'([A-Z])(\d+)=', r'\1\2\1', regex=True)

# '='가 포함되거나 공백인 행 제거
final_df = final_df[~final_df['Protein_Change'].str.contains('=', na=False) & final_df['Protein_Change'].str.strip().ne('')]
    
final_df = final_df[~final_df['Protein_Change'].str.contains(r'\?', na=False)]
# 'HGVSp_Short' 열에서 'fs' 뒤의 부분을 제거
final_df['Protein_Change'] = final_df['Protein_Change'].str.replace(r'fs.*', 'fs', regex=True)


# 변경된 결과 확인
print(final_df)

# 'SUBCLASS' 열을 제외한 train_df의 열들 가져오기
remaining_columns = train_df.drop(columns=['SUBCLASS','ID']).columns

# final_df의 'Hugo_Symbol'을 remaining_columns와 비교하여 필터링
final_df = final_df[final_df['Hugo_Symbol'].isin(remaining_columns)]

final_df

# 'HGVSp_Short' 열에서 'p.' 제거
final_df['Protein_Change'] = final_df['Protein_Change'].str.strip().str.replace('p.', '', regex=False)
# 그룹화하여 Protein_Change를 합치기
final_df = final_df.groupby(['Hugo_Symbol', 'Tumor_Sample_Barcode', 'SUBCLASS'], as_index=False).agg({
    'Protein_Change': lambda x: ' '.join(x.dropna().astype(str))  # NaN 값 제거 및 문자열로 변환 후 합치기
})
final_df
# 1. Hugo_Symbol 값의 고유 리스트 생성
hugo_symbols = train_df.drop(columns=['ID', 'SUBCLASS']).columns.tolist()

# 2. 새로운 DataFrame 초기화
result_rows = []  # 빈 리스트 초기화

# 3. Tumor_Sample_Barcode 및 SUBCLASS별로 데이터 처리
for (barcode, subclass), group in final_df.groupby(['Tumor_Sample_Barcode', 'SUBCLASS']):
    # 기본값을 가지는 행 생성 (모든 Hugo_Symbol은 'WT'로 설정)
    result_row = {hugo: 'WT' for hugo in hugo_symbols}
    
    # SUBCLASS와 ID 열 추가
    result_row['SUBCLASS'] = subclass
    result_row['ID'] = barcode  # Tumor_Sample_Barcode를 ID로 사용
    
    # 그룹 내의 각 Hugo_Symbol에 대해 Protein Change 값을 설정
    for _, row in group.iterrows():
        hugo_symbol = row['Hugo_Symbol']
        protein_change = row['Protein_Change']
        
        # 값이 이미 존재하는 경우, 띄어쓰기를 추가하여 값을 연결
        if result_row[hugo_symbol] != 'WT':
            result_row[hugo_symbol] += ' ' + protein_change  # 기존 값에 띄어쓰기를 추가하여 연결
        else:
            result_row[hugo_symbol] = protein_change  # 기본값을 설정
            
    result_rows.append(result_row)  # 결과 행을 리스트에 추가

# 4. 새로운 DataFrame 생성
result_df = pd.DataFrame(result_rows)

# 5. 빈 값은 'WT'로 대체
result_df.replace('', 'WT', inplace=True)

# 6. train_df의 열 추가 (없는 열을 'WT'로 채움)
for column in train_df.columns:
    if column not in result_df.columns:
        result_df[column] = 'WT'

result_df.to_csv("RESULT_FINAL.csv",index=False)
import pandas as pd

# 1. 데이터 로드 및 불필요한 열 제거
ALL_DATA_FINAL = pd.read_csv('ALL_DATA_FINAL.csv')
RESULT_FINAL = pd.read_csv('RESULT_FINAL.csv')
TRAIN = pd.read_csv('train.csv').drop(columns='ID')

# 2. TRAIN의 열 순서대로 ALL_DATA_FINAL과 RESULT_FINAL 정렬
ALL_DATA_FINAL_sorted = ALL_DATA_FINAL[TRAIN.columns]
RESULT_FINAL_sorted = RESULT_FINAL[TRAIN.columns]

# 3. 데이터프레임 병합 (기본적으로 행을 병합하려면 concat 사용)
FULL_DF = pd.concat([TRAIN, ALL_DATA_FINAL_sorted, RESULT_FINAL_sorted], ignore_index=True)


# 병합된 데이터 전처리
import pandas as pd
FULL_DF = pd.read_csv('FULL_DF.csv')
a_df = FULL_DF.drop(columns="SUBCLASS")
test_df = pd.read_csv('test.csv').drop(columns="ID")
# >가 들어간 데이터 분리
import pandas as pd
import re

def transform_variants(variant_str):
    # None 체크
    if pd.isna(variant_str):
        return variant_str  # None 또는 NaN인 경우 그대로 반환

    # 변이 문자열을 분리
    parts = variant_str.split(' ')
    transformed_parts = []
    
    # 패턴 정의
    pattern = r"(\d+)_(\d+)([A-Z*]+)>([A-Z*]+)"
    
    for part in parts:
        # 패턴에 매칭되는 경우
        match = re.match(pattern, part)
        if match:
            num1, num2, chars1, chars2 = match.groups()

            # 변이가 한 글자인 경우를 대비한 조건 처리
            result1 = chars1[0] + num1 + chars2[0] if len(chars1) > 0 and len(chars2) > 0 else ""
            result2 = ""
            
            if len(chars1) > 1 and len(chars2) > 1:
                result2 = chars1[1] + num2 + chars2[1]
            elif len(chars1) == 1 and len(chars2) == 1:
                result2 = chars1[0] + num2 + chars2[0]
            else:
                # 변이가 한 글자인 경우 처리
                result2 = chars1[0] + num2 + chars2[0]

            transformed = f"{result1} {result2}".strip()  # 빈칸이 남지 않도록 strip() 사용
            transformed_parts.append(transformed)
        else:
            # 패턴에 매칭되지 않는 부분은 그대로 유지
            transformed_parts.append(part)

    return ' '.join(transformed_parts)

# '>'가 포함된 변이를 변환
for col in a_df.columns:
    a_df[col] = a_df[col].apply(transform_variants)
    test_df[col] = test_df[col].apply(transform_variants)


# 멀티-핫-인코딩
#멀티핫 인코딩으로, trian_data에 대해 gene_mutation 구조로 열을 생성함 ex) A2M_R895R
import pandas as pd
from sklearn.preprocessing import MultiLabelBinarizer

# 변이 값을 리스트로 변환 (결측치는 빈 리스트로 대체)
def safe_split(x):
    if pd.isnull(x):
        return []
    return x.split()

# a_df에 대해 MultiLabelBinarizer를 사용하여 fit 수행
mlb_dict = {}  # 각 열에 대해 MultiLabelBinarizer를 저장하기 위한 딕셔너리
result_list_train = []  # 훈련 데이터의 결과를 저장할 리스트

for col in a_df.columns:
    multi_values = a_df[col].apply(safe_split)
    mlb = MultiLabelBinarizer()
    encoded = mlb.fit_transform(multi_values)
    
    # 열 이름을 col_변이 이름으로 설정
    col_encoded_df = pd.DataFrame(encoded, columns=[f"{col}_{variant}" for variant in mlb.classes_])
    
    # 데이터 타입을 boolean으로 변환
    col_encoded_df = col_encoded_df.astype(bool)
    
    result_list_train.append(col_encoded_df)
    
    # mlb를 딕셔너리에 저장
    mlb_dict[col] = mlb

# 결과 데이터프레임 합치기
multi_hot_df = pd.concat(result_list_train, axis=1)
# test_df에 대해 transform 수행
result_list_test = []  # 테스트 데이터의 결과를 저장할 리스트

for col in a_df.columns:
    multi_values_test = test_df[col].apply(safe_split)
    
    # 저장된 mlb를 사용하여 transform 수행
    encoded_test = mlb_dict[col].transform(multi_values_test)
    
    # 열 이름을 col_변이 이름으로 설정
    col_encoded_df_test = pd.DataFrame(encoded_test, columns=[f"{col}_{variant}" for variant in mlb_dict[col].classes_])
    
    # 데이터 타입을 boolean으로 변환
    col_encoded_df_test = col_encoded_df_test.astype(bool)
    
    result_list_test.append(col_encoded_df_test)

# 결과 데이터프레임 합치기
multi_hot_test_df = pd.concat(result_list_test, axis=1)
# 모델 학습
import pandas as pd
import lightgbm as lgb
from sklearn.preprocessing import LabelEncoder


label_encoder = LabelEncoder()
y_train = label_encoder.fit_transform(FULL_DF["SUBCLASS"])
lgb_train = lgb.Dataset(multi_hot_df, label=y_train)
# 모델 학습을 위한 파라미터 설정
params = {
    'objective': 'multiclass',
    'num_class': len(label_encoder.classes_),  # 클래스 수
    'metric': 'multi_logloss',
    'boosting_type': 'gbdt',
    'verbose': 1,
    'n_jobs': -1,
    'max_depth': -1,
    'num_leaves': 2000,
    'force_col_wise': True,
    'lambda_l1': 0.0,  # L1 정규화 제거
    'lambda_l2': 0.0 ,  # L2 정규화 제거
    'subsample' : 1,
}

# 모델 학습
model = lgb.train(params,
                  lgb_train,
                  num_boost_round=1000)
# 이후 feature importance를 이용하여 중요도가 높은 컬럼 제거를 반복 후 재학습
# 모델 예측
preds = model.predict(multi_hot_test_df)
# 예측된 클래스의 인덱스를 가져옴
predicted_classes = preds.argmax(axis=1)

# 인코딩된 클래스 인덱스를 원래 클래스 레이블로 변환
predicted_labels = label_encoder.inverse_transform(predicted_classes)
# 제출
submission = pd.read_csv('sample_submission.csv')
submission['SUBCLASS'] = predicted_labels
submission.to_csv('./FINALSUB.csv', encoding='UTF-8-sig', index=False)