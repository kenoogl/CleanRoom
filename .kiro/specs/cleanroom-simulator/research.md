# Research & Design Decisions

---
**Purpose**: クリーンルームシミュレータの技術設計を裏付ける調査結果と設計判断の記録
---

## Summary
- **Feature**: `cleanroom-simulator`
- **Discovery Scope**: New Feature（グリーンフィールド開発）
- **Key Findings**:
  - 参考実装（H2/src/）はモジュール分離型アーキテクチャを採用し、FLoops.jlで並列化を抽象化
  - Julia CFDエコシステムでは有限体積法+WENO再構成が標準的なアプローチ
  - CairoMakieが静的可視化出力に最適（PNG/SVG対応、ヘッドレス環境対応）

## Research Log

### 参考実装のアーキテクチャ分析

- **Context**: 既存のH2/src/熱解析コードからCFDソルバー設計のパターンを抽出
- **Sources Consulted**: `/Users/Daily/Development/H2/src/` の全ソースファイル
- **Findings**:
  - `Common.jl`: 共通定数、列挙型、WorkBuffers構造体、並列化ヘルパー
  - `BoundaryConditions.jl`: 境界条件の型定義と適用ロジック（structベース）
  - `NonUniform.jl`: 反復ソルバー（PBiCGSTAB, CG）と前処理
  - `RHS.jl`: 右辺項計算（境界条件のRHS寄与を分離）
  - `Zcoord.jl`: Z方向非等間隔格子座標の生成
  - 配列構造: (NX+2, NY+2, NZ+2)でゴーストセル込み
  - 並列化: `FLoops.jl`で`ThreadedEx()`/`SequentialEx()`を切り替え
- **Implications**:
  - モジュール分離パターンを継承
  - WorkBuffers構造体パターンを拡張して速度・圧力場を管理
  - FLoops.jlによる並列化抽象化を採用

### Julia CFDエコシステム調査

- **Context**: Julia CFD開発のベストプラクティスを確認
- **Sources Consulted**:
  - [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) - 有限体積法+WENO海洋シミュレータ
  - [XCALibre.jl](https://joss.theoj.org/papers/10.21105/joss.07441.pdf) - 非構造格子FVM
  - [CFD Julia](https://www.mdpi.com/2311-5521/4/3/159) - 教育用CFDモジュール
  - [CFD_Julia GitHub](https://github.com/surajp92/CFD_Julia) - WENO実装例
- **Findings**:
  - WENO-5よりWENO-3が低次元問題に適切（ステンシル幅2で十分）
  - Lax-Friedrichsフラックス分割は成分独立再構成に適合
  - Red-Black SORは並列化容易だがCG/BiCGSTABが収束性で優位
  - 非等間隔格子WENO係数は局所セル幅から動的計算
- **Implications**:
  - WENO3 + Lax-Friedrichs分割を採用（要件準拠）
  - Red-Black SOR必須、CG/BiCGSTABはオプション
  - 非等間隔格子係数はZcoord.jlパターンを参考に拡張

### 可視化ライブラリ調査

- **Context**: Julia内可視化機能のライブラリ選定
- **Sources Consulted**:
  - [Makie.jl公式ドキュメント](https://docs.makie.org/stable/)
  - [CairoMakieヒートマップ](https://docs.makie.org/dev/reference/plots/heatmap)
  - [Beautiful Makie例](https://beautiful.makie.org/dev/examples/2d/heatmaps/heatmap)
- **Findings**:
  - CairoMakie: 静的ベクター/ラスター出力（PNG, SVG, PDF）
  - GLMakie: GPU加速インタラクティブ表示
  - ヒートマップ: `heatmap(x, y, data; colormap=:viridis)`
  - コンター: `contourf` for filled contours
  - カラーバー: `Colorbar(fig, heatmap_plot)` で自動設定
  - 断面抽出: 3D配列のスライスを2D heatmapに直接渡せる
- **Implications**:
  - CairoMakieをデフォルト選択（ヘッドレス環境対応、高品質出力）
  - オプションでGLMakieインタラクティブモード

### SPH/バイナリ出力形式調査

- **Context**: V-Isio SPH形式とFortran unformatted互換性
- **Sources Consulted**: 要件ドキュメント内のSPHファイル仕様
- **Findings**:
  - リトルエンディアン、単精度Float32
  - Fortran unformatted形式: レコードマーカー（4バイトInt32）前後に配置
  - svType=1(スカラー), svType=2(ベクトル), dType=1(単精度)
  - 配列順序: column-major（Julia/Fortran互換）
- **Implications**:
  - Julia標準IOで直接書き込み可能
  - `write(io, htol(Int32(marker)))` でレコードマーカー

## Architecture Pattern Evaluation

| Option | Description | Strengths | Risks / Limitations | Notes |
|--------|-------------|-----------|---------------------|-------|
| Modular Monolith | 単一パッケージ内でモジュール分離 | シンプル、参考実装と整合 | 大規模化で複雑性増大 | **採用** |
| Hexagonal | ポート&アダプター抽象化 | テスト容易、外部依存分離 | Julia文化に馴染みにくい | 過剰設計 |
| Event-Driven | イベントソーシング | 状態追跡容易 | CFDには不適 | 却下 |

**選択理由**: 参考実装（H2/src/）のパターンを継承しつつ、CFDソルバー固有のモジュール（対流項、圧力ソルバー、境界条件等）を分離する「Modular Monolith」が最適。Juliaのモジュールシステムと`include` + `using`パターンに適合。

## Design Decisions

### Decision: モジュール構成

- **Context**: 12要件を満たすソルバーのコード構成
- **Alternatives Considered**:
  1. 単一ファイル実装 — シンプルだが保守困難
  2. モジュール分離（参考実装パターン） — 責務明確
  3. 外部パッケージ化 — 再利用性向上だが開発オーバーヘッド
- **Selected Approach**: モジュール分離（Option 2）
- **Rationale**: 参考実装で実証済み、チーム開発に適合
- **Trade-offs**: パッケージ化の再利用性は将来課題
- **Follow-up**: モジュール間依存の循環防止

### Decision: データ構造設計

- **Context**: 速度・圧力場と補助配列の管理
- **Alternatives Considered**:
  1. 個別配列（u, v, w, p 分離）
  2. 構造体（WorkBuffers拡張）
  3. NamedTuple
- **Selected Approach**: 構造体（CFDBuffers）
- **Rationale**: 参考実装のWorkBuffersパターンを継承、型安全性確保
- **Trade-offs**: 可変フィールドのミューテーション
- **Follow-up**: 配列サイズ検証をコンストラクタに実装

### Decision: 時間積分スキーム

- **Context**: Euler/RK2/RK4の切り替え
- **Alternatives Considered**:
  1. 条件分岐（if-else）
  2. 関数ディスパッチ（Val型）
  3. Strategy Pattern
- **Selected Approach**: Val型ディスパッチ
- **Rationale**: コンパイル時に分岐解決、ループ内オーバーヘッド回避
- **Trade-offs**: 型爆発の可能性
- **Follow-up**: ベンチマークで性能検証

### Decision: 可視化ライブラリ

- **Context**: 計算結果の画像出力
- **Alternatives Considered**:
  1. Plots.jl — 汎用だがバックエンド依存
  2. CairoMakie — 静的出力に最適
  3. PyPlot.jl — Python依存
- **Selected Approach**: CairoMakie + Package Extensions
- **Rationale**: Julia純正、PNG/SVG高品質出力、ヘッドレス対応
- **Trade-offs**: 初回コンパイル時間（Package Extensionsで軽減）
- **Follow-up**: `weakdeps`/`extensions`構成でCairoMakieを遅延読込

### Decision: 並列化バックエンド

- **Context**: 逐次/スレッド並列の切り替え
- **Alternatives Considered**:
  1. Threads.@threads直接使用
  2. FLoops.jl
  3. LoopVectorization.jl
- **Selected Approach**: FLoops.jl
- **Rationale**: 参考実装で実証済み、バックエンド抽象化
- **Trade-offs**: 依存パッケージ追加
- **Follow-up**: GPU対応（将来）はFolds.jl検討

## Risks & Mitigations

- **WENO3非等間隔格子の実装複雑性** — 要件に数式が明記済み、テスト駆動で検証
- **チェッカーボード不安定性** — セルフェイス内挿を要件で明確化済み
- **大規模格子でのメモリ不足** — dry_runモードでメモリチェック実装
- **SPHバイナリ互換性** — V-Isio実機で検証テスト必要
- **可視化の初回コンパイル遅延** — Package Extensions採用で非可視化ジョブの起動高速化

## References

- [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) — Julia有限体積法+WENO海洋シミュレータ
- [CFD Julia](https://www.mdpi.com/2311-5521/4/3/159) — 教育用CFDモジュール
- [Makie.jl Documentation](https://docs.makie.org/stable/) — Julia可視化エコシステム
- [FLoops.jl](https://github.com/JuliaFolds/FLoops.jl) — 並列ループ抽象化
- 参考実装: `/Users/Daily/Development/H2/src/` — 熱解析コード
