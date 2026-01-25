# 可視化ツール一覧（問題クラス別パターン）

内蔵可視化は汎用断面（速度・圧力）に限定し、問題特化の可視化はツールとして提供する。

## パターン一覧

### cavity
- **用途**: キャビティ流れの中心線プロファイル（Ghia形式）
- **入力**: 速度SPH（vel_#######.sph）
- **出力**: 中心線プロファイル図（`*_profile.png`）
- **ツール**: `tools/visualize_cavity.jl`

#### 例
```
julia tools/visualize_cavity.jl verification/cavity/output/vel_0001000.sph
julia tools/visualize_cavity.jl verification/cavity/output/vel 1000 2000
```

---

### backward_step
- **用途**: バックステップ流れの速度ベクトル＋圧力コンター
- **入力**: 速度SPH（vel_#######.sph）と圧力SPH（prs_#######.sph）
- **出力**: 速度ベクトル＋圧力コンター（`*_plot.png`）
- **ツール**: `tools/visualize_step.jl`

#### 例
```
julia tools/visualize_step.jl verification/backward_step/output/vel_0001000.sph
julia tools/visualize_step.jl verification/backward_step/output/vel 1000 2000
```

---

## 統一入口（dispatcher）

パターン指定で起動する統一入口を用意している。

```
julia tools/visualize.jl cavity verification/cavity/output/vel_0001000.sph
julia tools/visualize.jl cavity verification/cavity/output/vel 1000 2000
julia tools/visualize.jl backward_step verification/backward_step/output/vel_0001000.sph
julia tools/visualize.jl backward_step verification/backward_step/output/vel 1000 2000
```

---

## 追加パターンの方針
- 新しい問題クラスごとに専用ツールを追加する（例: `visualize_symmetric_channel.jl`）
- ツールはSPH出力を入力とし、可視化ロジックは自由に設計する
- 共通化が進んだ場合は統一入口（dispatcher）へパターン追加する
