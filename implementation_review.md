# 既存実装 精査レポート（仕様/設計準拠チェック）

## 対象
- 仕様書: `.kiro/specs/cleanroom-simulator/requirements.md`
- 設計書: `.kiro/specs/cleanroom-simulator/design.md`
- 実装: `src/` 以下一式

## まとめ（結論）
現行実装は**旧仕様（境界条件JSON・内部境界の扱い等）を前提にした実装が多く残っており、新仕様/設計に対して重要な不整合が複数あります。**
特に、境界条件スキーマ、圧力境界条件、内部境界の扱い、弱圧縮性項、平均化・可視化・SPHフォーマットなどで差分が顕著です。

---

## 重要度: 重大

### 1) 境界条件JSONのスキーマ／キーワード不一致（旧スキーマのまま）
- **仕様/設計**: `Opening`（flow_type=inlet/outlet）を含む新スキーマ。`Inflow/Outflow/Wall/...` を使用し、`dirichlet/neumann` は禁止。
- **実装**: `dirichlet/neumann/inlet` を受理、`inlets/outlets` 配列を前提。
- **影響**: 仕様準拠のJSONが読み込めず実行不可。
- 該当: `src/IO/InputReader.jl:259-377`、`src/BC/BoundaryConditions.jl:1-70`

### 2) internal_boundariesが「物体（mask=0）」として扱われていない
- **仕様/設計**: internal_boundaries は物体扱い（mask=0）、速度は物体内部速度。
- **実装**: 正規の境界条件として速度を上書きするのみで、mask=0化が無い。
- **影響**: 物体としての流束遮断・圧力処理が不一致。
- 該当: `src/IO/InputReader.jl:339-368`、`src/BC/BoundaryConditions.jl:780-815`、`src/BC/BoundaryConditions.jl:1142-1177`

### 3) 圧力境界条件が仕様と不一致（OutflowがDirichlet扱い）
- **仕様/設計**: Periodic以外は圧力Neumann固定。
- **実装**: `Outflow` を圧力Dirichletとして扱っている。
- **影響**: 圧力場の拘束条件が仕様と異なる。
- 該当: `src/BC/BoundaryConditions.jl:823-868`

### 4) 弱圧縮性連続式（M^2 dp/dt + ∂u_i/∂x_i = 0）が未実装
- **仕様/設計**: 弱圧縮性を考慮した連続式を前提。
- **実装**: その項を導入する箇所が存在しない（純粋な非圧縮 Fractional Step）。
- **影響**: 支配方程式レベルで仕様違反。
- 該当: `src/Physics/FractionalStep.jl`、`src/Physics/PressureSolver.jl`（全体）

---

## 重要度: 高

### 5) 前処理 sweep 数が仕様と不一致
- **仕様/設計**: Gauss–Seidel 4 sweep
- **実装**: `PRECONDITIONER_SWEEPS = 5`
- 該当: `src/Physics/PressureSolver.jl:14`

### 6) BiCGSTABが未実装
- **仕様/設計**: オプションとして選択時に適用。
- **実装**: 警告のみで実処理なし。
- 影響: `solver=BiCGSTAB`指定時に仕様違反。
- 該当: `src/Physics/PressureSolver.jl:430-446`

### 7) 圧力反復中の動的マスク対象がOutflowのみ
- **仕様/設計**: Outflow/Inflow/Opening を圧力反復中だけ mask=0。
- **実装**: `update_outflow_mask!` が Outflow のみ。
- 該当: `src/Physics/FractionalStep.jl:244-252`、`src/BC/BoundaryConditions.jl:1201-1237`

### 8) Outflow境界の粘性項ゼロが未対応
- **仕様/設計**: Outflowでは粘性項をゼロ。
- **実装**: DiffusionでOutflow専用の除外処理なし。
- 該当: `src/Physics/Diffusion.jl:1-200`

### 9) SPHフォーマットが仕様と不一致
- DZが非等間隔時に `dz[3]` を使用（仕様は Lz/Nz）
- ベクトル出力が**u,v,wのインターリーブ**（仕様は「u,v,wの順に配列を記述」＝ブロック列が想定）
- レコード構成（time/stepとorigin/pitch同レコード）
- 該当: `src/IO/SPHWriter.jl:61-129`

### 10) チェックポイントのエンディアン固定が未保証
- **仕様/設計**: Little Endian固定。
- **実装**: ネイティブエンディアンに依存。
- 該当: `src/IO/Checkpoint.jl:13-38`

---

## 重要度: 中

### 11) 時間平均の運用が未実装
- **仕様/設計**: Start_time_for_averaging以降、Welfordで平均更新＋SPH出力。
- **実装**: `update_time_average!` は定義のみ、呼び出し・平均出力が無し。
- 該当: `src/Fields.jl:153-169`、`src/CleanroomSolver.jl`（平均更新なし）

### 12) 可視化機能が仕様不足
- **仕様/設計**: velocity/pressure両方のコンター、vector、テキスト出力。
- **実装**: CairoMakie拡張は速度マグニチュード中心で、variables指定が無視される。
- 該当: `ext/CleanroomVisualizationExt.jl:35-150`、`src/IO/Visualization.jl`

### 13) 毎ステップの安定条件監視が未実装
- **仕様/設計**: CFL/拡散条件は毎ステップ監視。
- **実装**: 初期化時に拡散数のみチェック。`compute_cfl`は使われていない。
- 該当: `src/CleanroomSolver.jl:39-70`、`src/IO/Monitor.jl:88-126`

---

## 重要度: 低

### 14) 圧力平均引き戻しコメントが旧仕様表現
- コメントで「全境界Neumann」となっておりPeriodic例外が反映されていない。
- 該当: `src/Physics/PressureSolver.jl:62-66`

---

## 改修優先度（推奨）
1. **境界条件JSONスキーマの全面更新**（InputReader/BoundaryConditions/Monitorのcondition出力）
2. **internal_boundariesの物体化（mask=0）**と速度指定の整理
3. **圧力境界条件のNeumann統一（Periodic例外）**
4. **弱圧縮性項の導入方針の確定・実装**
5. **SPH/Checkpointフォーマットの仕様準拠**
6. **平均化・可視化・安定監視の実装**

---

## 付記
- 本レポートは「新仕様・設計が正」として現行実装を精査した結果です。
- 必要であれば、上記の差分を**具体修正パッチ**として提供できます。

