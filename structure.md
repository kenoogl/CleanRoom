# Structure Steering: cleanroom-simulator

## プロジェクト構造
Modular Monolithを採用し、責務ごとにモジュールを分離します。

```
src/
├── CleanroomSolver.jl      # メインエントリポイント・ドライバ
├── Common.jl               # 共通定数、型定義、無次元化ユーティリティ
├── Grid.jl                 # 格子生成、座標管理、Z方向非等間隔格子
├── Fields.jl               # 物理量配列（CFDBuffers）の管理
├── Physics/
│   ├── Convection.jl       # WENO3対流項離散化
│   ├── Diffusion.jl        # 2次精度中心差分拡散項
│   ├── Turbulence.jl       # LES (Smagorinsky) モデル
│   ├── PressureSolver.jl   # 圧力ポアソン方程式 (Red-Black SOR, CG, BiCGSTAB)
│   └── FractionalStep.jl   # 圧力-速度分離アルゴリズム
├── BC/
│   ├── BoundaryConditions.jl # 境界条件の定義と適用ロジック
│   └── Geometry.jl         # 物体形状定義、マスク生成
├── IO/
│   ├── InputReader.jl      # JSONパラメータ・境界条件読込
│   ├── SPHWriter.jl        # V-Isio SPHバイナリ出力
│   ├── Checkpoint.jl       # チェックポイント入出力
│   ├── Monitor.jl          # 計算モニター、history出力
│   └── Visualization.jl    # 可視化インターフェース（stub）
└── TimeIntegration.jl      # 時間積分スキーム (Euler, RK2, RK4)

ext/
└── CleanroomVisualizationExt.jl # CairoMakieによる可視化実装
```

## モジュール間依存関係
- `Common` はすべてのモジュールから参照される。
- `Physics` モジュール群は `Fields` と `Grid` に依存する。
- `FractionalStep` は `Physics` 内の各ソルバーを統合する。
- `CleanroomSolver` (Main) がすべてのモジュールをオーケストレーションする。

## 拡張設計
- **Package Extensions**: `CairoMakie` への依存は `ext/` に隔離し、可視化が不要な環境での起動高速化を図る。
- **HALO領域**: `Grid` および `Fields` は、将来的なMPI並列化を見据え、Z方向に2セルのHALO領域を確保する構造とする。