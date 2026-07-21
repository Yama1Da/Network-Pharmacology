import pandas as pd
import streamlit as st


@st.cache_data
def get_drug_targets(drug_name: str) -> list[str]:
    """读取 static/data/{drug_name}.csv，返回 'Gene' 列的靶点列表。"""
    df = pd.read_csv(f"static/data/{drug_name}.csv")
    return df["Gene"].tolist()


@st.cache_data
def get_disease_targets() -> list[str]:
    """返回硬编码的疾病靶点列表。"""
    return ["TP53", "AKT1", "TNF", "IL6", "VEGFA", "EGFR"]
