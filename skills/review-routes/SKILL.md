---
name: "review-routes"
description: "评审已有合成路线的拓扑、逐步反应先例、库存声明和证据缺口。用于比较 AiZynthFinder、PaRoutes 或结构化多步路线并生成专家复核队列。"
---

# 合成路线证据评审

## 能力

对已有合成路线执行确定性、非生成式证据评审：

- 读取 normalized route、AiZynthFinder JSON 或 PaRoutes v2 JSON；
- 校验单根、连通、无环的 `mol → reaction → mol` 交替路线树；
- 规范化 target、intermediate、terminal precursor 和每个单步反应；
- 用稳定 route/step hash 绑定 `curate-reactions` 与 `search-reactions` artifact；
- 传播 rejected、review、zero hit、timeout、source error 和 license 缺口；
- 检查库存快照、显式项目约束和重复路线；
- 分维度展示路线，不生成隐藏综合分数；
- 输出 weakest steps 和专家 review queue。

适用于：

- “检查这几条逆合成路线各有什么证据缺口”；
- “比较 AiZynthFinder 输出的几条路线”；
- “逐步检查路线中的反应先例”；
- “哪些路线步骤只有相似反应，没有精确先例”；
- “检查路线前体库存声明和来源许可”。

## 执行流程

1. 确认输入是已有路线 JSON，不接收 pickle，不调用路线生成模型。
2. 显式确认 `input_profile`、来源 SHA-256 和 `routes_fingerprint`。
3. 验证路线树并生成稳定 `route_signature/step_reaction_hash`。
4. 如提供第五、第六 Skill artifact，按 `route_id + step_id + step hash`
   验证后绑定，禁止按数组位置猜测。
5. 运行：

```bash
python scripts/review_routes.py \
  --input route-review-request.json \
  --output route-review.json
```

6. 校验：

```bash
python scripts/validate_output.py route-review.json
```

7. 报告路线 disposition、逐步 evidence level、库存/许可缺口、显式约束和
   review queue。

## 支持的输入 Profile

```text
normalized_route_v1
aizynthfinder_json
paroutes_v2_json
```

首版最多 20 条路线、每条最多 50 步、总节点最多 5000。

## 强制边界

- 固定 `rdkit==2025.9.2`；
- 核心不依赖 AiZynthFinder、Syntheseus、OpenAI4S 或模型权重；
- 只读取 JSON，禁止 pickle 或可执行反序列化；
- 不生成、补全或自动修复路线；
- 不运行 forward/round-trip、可行性、条件、yield 或安全模型；
- backend rank/score 原样保留，不跨 backend 归一化；
- exact、transformation、similar、component、zero hit 和 provider failure
  必须分开；
- 条件和产率仅为来源报告证据，不自动迁移到目标步骤；
- 库存必须关联快照；route export 的 `in_stock` 只标记为来源报告；
- 缺失许可保持 `null`，不得猜测；
- 不删除重复路线，只分组；
- 默认只输出 `dimensions_only`，不计算综合总分；
- `ready_for_expert_review` 仅表示证据包可交给专家，不表示路线可行、安全、
  最优、可放大或可直接实验；
- 不访问网络，不调用商业库存、路线生成器或远程反应数据库。

## 与其他 Skill 的关系

```text
curate-reactions
单步反应结构、角色、产率和质量状态
        ↓
search-reactions
单步先例、相似度、条件、产率、来源和许可
        ↓
review-routes
多步拓扑、逐步证据覆盖、最弱步骤和专家复核队列
```

完整输入输出、状态、规则和科学边界见
`references/输入输出与科学边界.md`。
