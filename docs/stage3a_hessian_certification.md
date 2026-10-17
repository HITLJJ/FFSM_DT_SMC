# Stage 3A：Taylor–Hessian 严格误差界认证

## 结论

Stage 3A 已完成。在与 Stage 2 完全相同的 6 个固定状态、3 个参数盒和 2 个控制通道组成的 36 个案例中，完整参考力矩映射

\[
F(z,\theta)=M_r(q,\theta)\ddot q_r+C_r(q,\dot q,\theta)\dot q_r
\]

均获得 `CERTIFIED` Taylor–Hessian 界，认证失败数为 0。四个额外压力状态也全部完成认证。

新界相对 Stage 2 严格基线的缩小倍数中位数为 66.6593，最大值为 704.6545；但相对经验搜索最大误差仍有 4.4876e4 的中位保守倍数和 2.1119e6 的最大保守倍数。因此本阶段解决了“Taylor Hessian 是否可以严格认证”的问题，但没有达到接近经验最坏误差的低保守目标。

## 版本与边界

- 起点提交：`58963429a10a3b3216cb6704601a5e4c7f134607`。
- 分支：`stage3-hessian`。
- MATLAB：R2025b（25.2），Symbolic Math Toolbox 25.2。
- 没有安装或使用 INTLAB、`intval`、`infsup` 等外部区间库。
- 没有修改 `reference_torque_map.m` 或 `src/dynamics/` 下的 Stage 1 核心动力学。
- 没有实现 SMC、DOB、在线辨识或数字孪生更新。

认证范围是：对每个给定固定状态 \(z\)，覆盖完整的正惯性参数盒 \(\Theta\)。状态、几何和参数盒均显式记录。经验参数点只用于评估，认证器不读取所谓“真实被控对象参数”。

## Hessian 推导路径

记

\[
K=A^{-1},\qquad S=B_c^TKB_c,\qquad M_r=D-S.
\]

由于完整质量矩阵对六个惯性参数仿射，\(A_{jk}=B_{c,jk}=D_{jk}=0\)，并有

\[
K_j=-KA_jK,
\]

\[
K_{jk}=KA_kKA_jK+KA_jKA_kK.
\]

一阶约化惯性参数导数为

\[
M_{r,j}=D_j-B_j^TKB-B^TKB_j+B^TKA_jKB.
\]

二阶导数由八个 Schur 乘积词组成：

\[
\begin{aligned}
M_{r,jk}={}&-B_j^TKB_k-B_k^TKB_j
+B_j^TKA_kKB+B_k^TKA_jKB\\
&+B^TKA_jKB_k+B^TKA_kKB_j\\
&-B^TKA_kKA_jKB-B^TKA_jKA_kKB.
\end{aligned}
\]

认证器对上述每个乘积词的每个因子应用关节坐标乘积法则，并使用

\[
K_{,q_\ell}=-KA_{,q_\ell}K
\]

生成 \(M_{r,j,q_\ell}\) 和 \(M_{r,jk,q_\ell}\) 的完整范数覆盖。随后按 Stage 1 相同的 Christoffel 指标约定构造 \(C_{r,j}\) 与 \(C_{r,jk}\)，最终得到

\[
F_{,j}=M_{r,j}\ddot q_r+C_{r,j}\dot q_r,
\qquad
F_{,jk}=M_{r,jk}\ddot q_r+C_{r,jk}\dot q_r.
\]

实施前的逐项数学—代码对应关系见 `docs/stage3a_math_code_map.md`。

## 严格数值覆盖方法

### 固定构型三角函数包络

`fixed_state_block_norms.m` 将 binary64 关节角转换为精确二进制有理数。`sin` 和 `cos` 使用 18 项精确有理 Taylor 多项式，并分别附加 Lagrange 余项

\[
|R_{\sin}|\le \frac{|x|^{36}}{36!},\qquad
|R_{\cos}|\le \frac{|x|^{35}}{35!}.
\]

质量 Jacobian、关节坐标导数以及 \(A_j,B_j,D_j,A_{j,q},B_{j,q},D_{j,q}\) 均通过精确有理端点矩阵运算传播。普通 double 端点没有被当成区间证书。

### 逆矩阵与 Schur 乘积

正惯性项给出

\[
A\succeq\operatorname{diag}(m_0,m_0,J_0+J_1+J_2),
\]

因此

\[
\|A^{-1}\|_2\le
\frac{1}{\min\{\underline m_0,
\underline J_0+\underline J_1+\underline J_2\}}.
\]

每个 Schur 乘积词均使用精确有理范数系数、次乘性和乘积法则独立传播。最终 double 输出由符号比较验证其不小于精确表达式，必要时向上移动 ULP。

该实现不需要参数盒划分，`subbox_count=1`。因此不存在被丢弃的失败子盒。若输入盒无效、Symbolic Toolbox 不可用或任一认证步骤抛出异常，函数返回 `FAILED` 和具体 `failure_reason`。

## 线性项和 Hessian 余项

为了保持完整计算链严格，线性项没有使用有限差分 Jacobian 点值，而是使用名义参数处解析梯度的严格分量上界：

\[
\beta_{i,\mathrm{lin}}
\le \sum_j \overline{|F_{i,j}(\hat\theta)|}\,r_j.
\]

对每个通道构造全部 36 个 Hessian 元素包络

\[
|H_{i,jk}(\theta)|\le U_{i,jk},\quad \forall\theta\in\Theta.
\]

比较三种严格余项：

1. 逐元素：
   \[
   R_i^{elem}=\frac12\sum_{j,k}U_{i,jk}r_jr_k.
   \]
2. 原始参数谱/Frobenius：
   \[
   R_i^{raw}=\frac12\|U_i\|_F\|r\|_2^2.
   \]
3. 无量纲归一化谱/Frobenius。

令 \(S=\operatorname{diag}(s_1,\ldots,s_6)\)，\(s_j=\hat\theta_j\)，以及

\[
\xi=S^{-1}(\theta-\hat\theta).
\]

则

\[
J_\xi=J_\theta S,
\qquad
H_\xi=S^TH_\theta S.
\]

归一化只改变坐标表达，不改变参数盒或力矩单位。逐元素余项在相同对角缩放下保持不变；统一谱/Frobenius 界则因参数尺度更均衡而收紧。

## 解析导数独立交叉验证

`reference_torque_parameter_derivatives.m` 按上述 Schur 公式直接计算一阶、二阶参数导数，并以复步长计算其关节坐标导数。它与 Stage 2 的中心有限差分 `parameter_jacobian.m`、`parameter_hessian.m` 不共享 Hessian 展开。

名义弯曲状态下：

- Jacobian Frobenius 差：2.2730e-11；
- Hessian Frobenius 差：9.3166e-9；
- Hessian 参数指标对称误差：9.8531e-17。

有限差分和复步长仅用于实现交叉检查，不参与认证常数生成。

## 主实验结果

| 指标 | 结果 |
|---|---:|
| 主通道案例 | 36 |
| `CERTIFIED` | 36/36 |
| 认证失败 | 0 |
| Stage 3A 界范围 | 50.4908–576431.6841 N·m |
| Stage 2/Stage 3A 缩小倍数，中位/最大 | 66.6593 / 704.6545 |
| Stage 3A/经验最大值，中位/最大 | 44876.2773 / 2111921.43 |
| Hessian 余项占总界，中位/最大 | 73.7025% / 93.3910% |
| Hessian 余项占总界，最小 | 35.9488% |
| 单次盒认证部分最大耗时 | 2.7781 s |
| 固定状态系数认证最大耗时 | 29.1459 s |

所有经验最大值均来自角点、固定种子随机点和多起点 `fmincon`，只用于反例搜索。36 个主案例没有发现严格界违反。

### 余项方法比较

| 方法 | 完整界中位数 [N·m] | 结论 |
|---|---:|---|
| 逐元素 | 5004.5134 | 36/36 案例最紧 |
| 原始参数统一谱界 | 14988.0297 | 最宽 |
| 归一化参数统一谱界 | 10350.1626 | 比原始尺度谱界缩小约 30.9%，但仍宽于逐元素界 |

这是一项负面与正面并存的结果：归一化确实改善了统一谱界的尺度条件，但没有超过逐元素方法。

## 压力测试

全部压力案例使用 ±20% 参数盒并保留在结果中：

| 压力项 | 状态 | 通道 1 界 [N·m] | 通道 2 界 [N·m] |
|---|---|---:|---:|
| 参数盒边界反例搜索 | `CERTIFIED` | 107376.8188 | 104910.0150 |
| 接近折叠的运动学构型 | `CERTIFIED` | 810.3381 | 680.9081 |
| 大参考速度 | `CERTIFIED` | 1771187.0006 | 1673696.1385 |
| 大参考加速度 | `CERTIFIED` | 60309.1347 | 60657.5593 |

“接近奇异”仅指运动学折叠压力状态；正质量基座保证本研究中的 \(A\) 仍具有所用的正定下界。大参考速度案例最不利，反映 Christoffel/Hessian 范数传播仍是主要瓶颈。

## 测试与复现

- Stage 1：10/10 通过；独立符号一致性检查通过。
- Stage 2：10/10 通过。
- Stage 3A：5/5 通过。
- 没有修改既有测试阈值。

运行：

```matlab
cd('D:\LJJ\FFSM_DT_SMC');
addpath(genpath(pwd));
run_stage3a_tests;
summary = run_stage3a_hessian_certification;
```

结果目录 `results/stage3a/` 包含：

- `stage3a_summary.csv`：36 个主通道案例；
- `stress_summary.csv`：全部压力案例；
- `stage3a_raw.mat`：配置、原始记录和失败列表；
- `summary.txt`；
- 6 组 PNG 与可编辑 FIG 图。

## Q1–Q6

**Q1：是否完成完整参考力矩映射的全域 Hessian 认证？**  是。对每个固定 \(z\)，全部 \(F_i\) 的 6×6 参数 Hessian 元素均在完整参数盒上获得严格包络；Schur 补、关节坐标导数和 Coriolis 二阶参数导数均包含在内。

**Q2：是否获得严格 Taylor 误差界，而非候选界？**  是。36 个主案例和 4 个压力案例均为 `CERTIFIED`。Stage 2 的采样 Hessian 方法仍保留为 `CANDIDATE`，没有被重新命名。

**Q3：相对于 Stage 2 严格基线，保守性降低了多少？**  上界宽度缩小倍数的中位数为 66.6593，最大值为 704.6545；最差案例仍缩小约 8.0709 倍。但相对经验误差仍高度保守，不能宣称已经得到低保守界。

**Q4：哪一种余项算法最紧？它的严格性如何保证？**  逐元素算法在 36/36 案例中最紧。其严格性来自每个 \(U_{i,jk}\) 的精确有理三角包络、Schur 乘积词范数传播、Christoffel 传播和受控向上舍入；不是来自采样最大值。

**Q5：如果未达到低保守目标，主要瓶颈是什么？**  主要瓶颈是区间/范数依赖性，其次是逆矩阵统一下界和 Christoffel 中多次绝对值求和。参数尺度确实影响统一谱界，但归一化后仍不如逐元素方法，因此参数尺度不是当前首要瓶颈。余项占总界中位 73.70%，说明 Hessian 包络仍是核心收紧对象。

**Q6：是否具备进入 Stage 3B 的条件？**  有条件具备。严格 Hessian 和完整 Taylor 界已经建立，可进入 Stage 3B 研究子盒自适应划分、经验证逆矩阵包络及依赖性削减；但在进一步收紧前，不建议把当前数值直接用于 SMC 鲁棒增益，也不应进入闭环控制实现。

Stage 3A 到此停止。
