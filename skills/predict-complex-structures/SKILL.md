---
name: "predict-complex-structures"
description: "基于科研证据、已审查蛋白质与配体输入及 Protenix 生成并校验生物分子复合物结构预测。用于蛋白质单体、多链或蛋白质–小分子复合物的可追溯预测与专家复核准备。"
---

# 生物分子复合物结构预测

把经过来源核对和化学结构审查的生物分子输入转换为官方兼容的
Protenix 任务，在独立 NVIDIA GPU Runtime 中执行预测，并生成可校验、
可追溯的专家复核包。

本 Skill 生成的是结构假设，不是实验结构、结合证明、活性结论或治疗建议。

## 执行流程

### 1. 发现科研证据

先由宿主 Agent 调用已经配置的 `dataPro_search` MCP Tool，检索靶蛋白、
配体、实验结构、结合位点、突变和相关研究。

- Skill 内不得实现 DataPro Client；
- 不读取、要求、打印或保存 Agent Plan Key；
- 将 DataPro 原始结果保存为 JSON Artifact，再交给本 Skill；
- 题名、URL、摘要、作者、日期和 DOI 原样保留；
- 搜索命中和摘要只属于候选证据，不能写成全文已核验；
- DataPro 不可用时不伪造结果，后续输入准备必须写
  `evidence_status=not_run`。

### 2. 审查蛋白质和配体输入

蛋白质链必须提供明确序列、链 ID、物种和来源。UniProt、RCSB PDB 等
主数据由宿主 Agent 或用户提供，本 Skill 不使用 DataPro 摘要替代序列。

配体按输入类型复用已有化学 Skill：

1. 名称、CAS RN、PubChem CID、ChEMBL ID 等先调用
   `resolve-chemical-identities`；
2. SMILES、InChI、SDF 或已解析结构再调用
   `standardize-chemical-structures`；
3. 仅消费两个 Skill 保存的 JSON Artifact，不复制其算法；
4. `review_required` 继续传播，但用户确认具体计算形式后可以生成模型输入；
5. `rejected`、`not_run` 或 `error` 必须阻止 Protenix；
6. `parent` 只保留为派生比较视图，绝不静默代替实际配体；
7. `compute-molecular-features` 和
   `search-and-curate-chemical-libraries` 仅在用户明确要求候选筛选时调用。

### 3. 准备 Protenix 输入

在本 Skill 的轻量隔离环境中安装固定输入/验证依赖：

```bash
python -m pip install -r scripts/requirements.txt
```

准备结果目录：

```bash
python scripts/prepare_input.py \
  --request request.json \
  --result-dir complex-result
```

该命令生成：

- `input_manifest.json`；
- `evidence.json`；
- `ligand_handoff.json`；
- `protenix_input.json`；
- 初始 `run_manifest.json`。

`protenix_input.json` 顶层严格为官方要求的任务数组，内部只使用
`proteinChain`、`dnaSequence`、`rnaSequence`、`ligand`、`ion`、
`covalent_bonds` 和 `constraint` 等官方字段。文件型输入会复制到结果目录的
`inputs/`，模型输入只记录相对路径。请求引用的 Artifact、MSA、模板、配体
结构和参考结构必须位于 `request.json` 所在目录内；文件型科学输入必须提供
SHA-256，禁止绝对路径越界或 `..` 逃逸该目录。

### 4. 在独立 Runtime 执行 Protenix

Protenix 及其 CUDA、模型权重、MSA/模板工具与现有化学 Skill 使用独立
Runtime。不得修改任一侧 RDKit 版本来强行合并环境。
首版只接受 reference 中列出的受控模型；需要额外 ESM/ISM 权重的模型不自动
启用。

GPU Runtime 必须先生成 `runtime-assets.json`，记录
`PROTENIX_ROOT_DIR`、实际 Protenix 版本、模型 checkpoint 及必需缓存的
仓库相对路径和 SHA-256。默认禁止 Protenix 隐式下载。

在已配置好权重和缓存的 GPU Runtime 生成清单：

```bash
python scripts/build_runtime_assets.py \
  --protenix-root-dir /managed/runtime/protenix-data \
  --model-name protenix_base_default_v1.0.0 \
  --output runtime-assets.json
```

GPU Runtime 中运行：

```bash
python scripts/run_protenix.py \
  --result-dir complex-result \
  --protenix-command protenix \
  --runtime-assets-manifest runtime-assets.json \
  --timeout-seconds 3600
```

每次执行必须使用刚由 `prepare_input.py` 创建、从未运行过的新结果目录；失败、
取消或超时后若要重试，重新 prepare 到另一目录。Runner 拒绝含旧日志、旧报告
或旧预测文件的目录，避免把陈旧结果误认成本次输出。

只有用户明确确认下载范围、磁盘和网络影响时，才可用
`--allow-download` 代替资产清单。

Runner 实际执行的命令遵循官方 CLI：

```bash
protenix pred \
  -i protenix_input.json \
  -o protenix_output \
  -s 101 \
  -n protenix-v2 \
  --use_msa false \
  --use_template false \
  --use_rna_msa false
```

实际参数来自 `run_manifest.json`，不得用上面的示例覆盖用户任务。

Runner 必须：

- 检查 Protenix CLI 和 NVIDIA GPU；
- 默认验证 checkpoint、运行缓存、Protenix 版本和文件 SHA-256；
- 未显式 `--allow-download` 时禁止启动可能下载资产的官方 CLI；
- 启动子进程前移除 Key、Token、Authorization、Cookie、Secret 和 Password
  类环境变量；
- 保存经过脱敏的 stdout、stderr 和 exit code；
- 支持超时、取消文件和进程组终止；
- 保留官方生成的 CIF/mmCIF 与 summary/full confidence JSON；
- 将官方 `_wounresol.cif` 作为辅助结构单独登记，不与主样本重复配对；
- 不通过日志文本猜测成功；
- 必需文件缺失或 Validator 失败时不得写 `completed`。

缺少 GPU、CLI、权重或依赖时输出 `not_run`，并说明缺少的资源。不得生成
伪 CIF 或伪 confidence。

### 5. 独立校验结果

无论是否执行 GPU，都运行：

```bash
python scripts/validate_output.py complex-result
```

Validator 检查：

- 五个基础 JSON Artifact 的结构和确定性指纹；
- CIF/mmCIF 可解析性、链 ID、实体和原子数；
- 已确认配体的重原子数与元素组成；
- summary confidence 的有限数值和复合物链对字段；
- 官方可选的逐链/链对置信度、碰撞标记和 recycle 次数；
- 结构与 confidence 的 job、seed 和 sample 必须一一对应；
- 上游状态是否完整传播；
- 输出中是否包含 Secret 或主机敏感路径；
- `completed` 是否同时满足 exit code 0、必需文件存在和全部合同检查。

当前 Validator 不独立证明共价键或软约束已经被预测结果满足；只要输入含
`covalent_bonds` 或 `constraint`，结果必须保持 `review_required`，由领域专家
结合原子距离、化学合理性和实验背景复核。

Validator 生成：

- `validation_report.json`；
- `report.html`。

## 状态

运行状态：

```text
completed
partial
not_run
error
```

审查结论：

```text
ready_for_review
review_required
rejected
```

`ready_for_review` 只表示结果包可以交给结构生物学、计算化学或药物研发
人员复核，不表示预测正确或可以直接实验。

## 强制科学边界

- 预测结构不是实验结构；
- ranking score、pLDDT、pTM、ipTM 和 chain-pair ipTM 不是 Kd、IC50、
  结合自由能、药效或实验准确率；
- 高置信度不能排除错误构象、错误质子化、错误配体形式或错误生物学组装；
- 多 seed 一致只表示当前模型与参数下较稳定，不证明真实结构正确；
- pocket/contact 是软约束；满足约束不能作为独立实验验证；
- DataPro 命中或摘要不等于论文全文已经核验；
- 不自动声称“已证明结合”“已确认活性”“可直接用于实验”或“具有治疗效果”；
- 不生成实验执行、安全、临床或治疗结论；
- 模型结果必须由相关领域专家复核。

完整合同见[输入输出与科学边界](references/输入输出与科学边界.md)。
