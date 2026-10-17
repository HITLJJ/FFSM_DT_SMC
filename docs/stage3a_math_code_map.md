# Stage 3A 数学项与代码路径对应表（实施前）

起点提交：`58963429a10a3b3216cb6704601a5e4c7f134607`；工作分支：`stage3-hessian`。

本表在 Stage 3A 算法代码修改前建立。记

\[
K=A^{-1},\qquad S=B_c^T K B_c,\qquad M_r=D-S.
\]

完整质量矩阵对六个惯性参数仿射，因此
\(A_{jk}=B_{c,jk}=D_{jk}=0\)。几何参数和控制相关状态 \(z\) 固定。

| 数学对象 | 展开或覆盖路径 | 计划代码责任 |
|---|---|---|
| \(A_j,B_j,D_j\) | 完整质量矩阵仿射参数系数块 | `affine_mass_block_derivatives.m` 提供数值解析系数，用于交叉验证；认证器使用 `analytic_block_norms.m` 的精确有理全局系数界 |
| \(K_j\) | \(-KA_jK\) | `reference_torque_parameter_derivatives.m` 直接矩阵计算；认证器用 \(\|K\|^2\|A_j\|\) |
| \(K_{jk}\) | \(KA_kKA_jK+KA_jKA_kK\) | 同上；认证器逐乘积项传播，不把采样值用于覆盖 |
| \(M_{r,j}\) | \(D_j-B_j^TKB-B^TKB_j+B^TKA_jKB\) | 解析导数函数计算；认证器对四项分别作精确范数上界 |
| \(M_{r,jk}\) | \(-B_j^TKB_k-B_k^TKB_j+B_j^TKA_kKB+B_k^TKA_jKB+B^TKA_jKB_k+B^TKA_kKB_j-B^TKA_kKA_jKB-B^TKA_jKA_kKB\) | 解析导数函数计算实际矩阵；`certified_hessian_bound.m` 将八个乘积词逐项覆盖 |
| \(\partial M_{r,j}/\partial q_\ell\) | 对一阶 Schur 乘积词使用乘积法则；\(K_{,\ell}=-KA_{,\ell}K\) | 解析值以独立复步长关节导数生成；认证器把每个词的每个可微因子展开并求和 |
| \(\partial M_{r,jk}/\partial q_\ell\) | 对上述八个 Hessian 乘积词逐因子使用乘积法则 | 同上；认证路径使用 \(A_{j,q_\ell},B_{j,q_\ell}\) 和盒内 \(A_{,q_\ell},B_{,q_\ell}\) 的精确全局界 |
| \(C_{r,j}\) | Christoffel 公式作用于 \(M_{r,j,q_\ell}\) | `reference_torque_parameter_derivatives.m` 按 Stage 1 相同指标约定组装 |
| \(C_{r,jk}\) | Christoffel 公式作用于 \(M_{r,jk,q_\ell}\) | 同上；认证器由导数矩阵谱范数界得到每个 Coriolis 元素界 |
| \(F_{,j}\) | \(M_{r,j}\ddot q_r+C_{r,j}\dot q_r\) | 解析导数函数返回实际 Jacobian；认证器返回其逐参数绝对值上界 |
| \(F_{,jk}\) | \(M_{r,jk}\ddot q_r+C_{r,jk}\dot q_r\) | 解析导数函数返回实际 Hessian；认证器返回全部 \(U_{i,jk}\) |
| 线性支持项 | \(\sum_j |F_{i,j}(\hat\theta)|r_j\) | 若实际 Jacobian 没有经验证的点值区间，则使用严格的 \(|F_{i,j}|\) 解析上界，绝不以有限差分值标记认证 |
| 逐元素余项 | \(\frac12\sum_{j,k}U_{i,jk}r_jr_k\) | `certified_hessian_bound.m` 的首选严格余项 |
| 统一谱余项 | \(\frac12L_i\|r\|_2^2\)，\(L_i\le\|U_i\|_F\) | 同一认证器计算并与逐元素方法比较 |
| 归一化参数 | \(\xi=S^{-1}(\theta-\hat\theta)\)，\(H_\xi=S^TH_\theta S\) | 比较原尺度与 \(s_j=\hat\theta_j\) 的严格余项；不改变物理参数盒或 N·m 单位 |

独立交叉验证路径为：解析 Schur 参数导数与 `parameter_hessian.m` 的中心有限差分结果比较；两者不共享 Hessian 展开。有限差分只用于发现实现错误，不参与 `CERTIFIED` 常数生成。
