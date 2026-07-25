import pandas as pd
import streamlit as st

# 读取药物靶点列表（测试假数据）
@st.cache_data
def get_drug_targets(drug_name: str) -> list[str]:
    """读取 static/data/{drug_name}.csv，返回 'Gene' 列的靶点列表。"""
    df = pd.read_csv(f"Static/data/{drug_name}.csv")
    return df["Gene"].tolist()

# 读取疾病靶点列表（测试假数据）
@st.cache_data
def get_disease_targets(disease_name: str) -> list[str]:
    """读取 static/data/{disease_name}.csv，返回 'Gene' 列的靶点列表。"""
    df = pd.read_csv(f"Static/data/{disease_name}.csv")
    return df["Gene"].tolist()
