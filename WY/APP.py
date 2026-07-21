import streamlit as st
import UTILS

st.set_page_config(page_title="网络药理学平台", layout="centered")
st.title("🧬 网络药理学分析平台")

# ---------- 侧边栏 ----------
with st.sidebar:
    st.header("分析参数")
    drug_name = st.text_input("药物名称", value="aspirin",
                               help="输入小分子药物名称（如 aspirin）")
    analyze_btn = st.button("开始分析", type="primary")

# ---------- 主区域 ----------
if analyze_btn:
    with st.spinner(f"正在查询 {drug_name} 的靶点信息…"):
        drug_targets = UTILS.get_drug_targets(drug_name)

    st.success(f"✅ 分析完成")
    st.write(f"药物 **{drug_name}** 的靶点数量：**{len(drug_targets)}**")

    with st.expander("查看靶点列表"):
        st.write(drug_targets)
else:
    st.info("👈 请在侧边栏输入药物名称，然后点击「开始分析」")
