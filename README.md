# SciAgent Chemistry Skills

[![CI](https://github.com/3494036618-eng/sciagent-chemistry-skills/actions/workflows/ci.yml/badge.svg)](https://github.com/3494036618-eng/sciagent-chemistry-skills/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

一组面向科研 Agent 的可审计化学 Skill，覆盖小分子身份、结构、二维特征、分子库治理、反应整理、反应检索和已有合成路线评审。

> English summary: Auditable chemistry skills for scientific agents, with deterministic computation, provenance, controlled status propagation, and explicit scientific limits.

## 状态

- 版本：`0.1.0-alpha.1`
- 阶段：工程候选版
- 公开测试：306 项
- 化学专家验收：尚未完成
- 真实用户验收：尚未完成
- 生产可用声明：无

测试通过只证明已覆盖的代码合同在指定环境中可复现，不证明化学结论、实验可行性、安全性或业务适用性。

## Skills

| Skill | 作用 | 主要边界 |
|---|---|---|
| `resolve-chemical-identities` | 保守解析名称和标识符，保留多源候选、歧义与冲突 | 在线名称会发送给外部服务；不接入商业 CAS 数据库 |
| `standardize-chemical-structures` | 解析、标准化、提取 parent 并执行结构 QC | parent 不是物理样品；不确认身份或活性 |
| `compute-molecular-features` | 计算受控二维描述符和 Morgan、RDKit、MACCS 指纹 | 描述符不是实验值；相似指纹不代表活性相同 |
| `search-and-curate-chemical-libraries` | 相似性、子结构、聚类、多样性选择和只读治理 | 不自动删除、合并或改写原始化合物库 |
| `curate-reactions` | 非破坏性整理 reaction SMILES、ORD 和反应表格 | 守恒与映射仅作诊断，不证明反应正确 |
| `search-reactions` | 按 ID、组分、SMARTS 或整体反应相似度检索先例 | 先例不等于路线可行或条件推荐 |
| `review-routes` | 审查已有路线拓扑、逐步证据、库存声明和证据缺口 | 不生成逆合成路线，不判断路线安全、最优或可实验 |

典型衔接：

```text
resolve identities
→ standardize structures
→ compute features
→ search/curate libraries

curate reactions
→ search reaction precedents
→ review existing routes
```

## 设计原则

1. 原始记录、派生记录和失败记录分别保留。
2. 确定性化学计算交给固定版本的专业工具。
3. 工具、版本、参数、来源、状态和结果指纹进入输出合同。
4. `rejected`、`review_required`、零命中、超时和来源错误不得静默丢失。
5. 结构相似、数据库命中和反应先例不得自动升级为科学结论。
6. Skill 不依赖 SciAgent 平台，可由其他 Agent runtime 独立调用。

每个 Skill 的完整输入输出合同和科学边界位于其 `references/` 目录。

## 环境

- Python `3.11` 或 `3.12`
- macOS 或 Linux
- RDKit `2025.9.2`
- ChEMBL Structure Pipeline `1.2.4`
- ORD Schema `0.8.3`

建议创建隔离环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

运行仓库结构和公开边界检查：

```bash
python scripts/validate_repository.py
```

各 Skill 的触发语义、CLI 参数和输出状态见对应 `SKILL.md`。

## 网络与数据

大部分 Skill 离线运行。以下能力可访问外部服务：

- `resolve-chemical-identities`：OPSIN、PubChem、ChEMBL、UniChem；
- `search-reactions`：Open Reaction Database 官方 API。

敏感、保密或未公开的名称和结构不应发送给第三方服务。仓库不包含 ORD、PaRoutes、ChEMBL、PubChem 或 Agent Plan 的真实数据导出和评测结果。运行时取得的数据仍受原始来源条款约束。

第三方软件、服务和数据边界见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。

## 贡献与安全

- 贡献规则：[CONTRIBUTING.md](CONTRIBUTING.md)
- 安全报告：[SECURITY.md](SECURITY.md)
- 行为准则：[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

涉及科学规则、阈值、状态升级或来源解释的修改必须附一手依据和回归测试。

## 许可证

项目原创代码和文档以 [Apache License 2.0](LICENSE) 发布，Copyright 2026 yutong。

第三方依赖和运行时数据不因此转为 Apache-2.0；其许可证和归属要求分别适用。
