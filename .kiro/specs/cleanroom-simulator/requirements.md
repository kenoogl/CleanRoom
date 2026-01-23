# Requirements Document

## Introduction
クリーンルーム内微小粒子の発生源特定を目的とした3次元CFDソルバー。弱圧縮性を考慮した非圧縮性Navier-Stokes方程式をLESモデルで解き、後段でBackward FTLE + Backward LPT解析を行うための速度場データを生成する。開発言語はJulia。

---

## Functional Requirements

### Requirement 1: 計算格子生成
**Objective:** As a シミュレーションエンジニア, I want 矩形計算領域に直交格子を生成できる, so that 3次元CFD計算の空間離散化ができる

#### Acceptance Criteria
1. The Solver shall 矩形計算領域に直交格子を生成する
2. The Solver shall XY方向に等間隔格子（Lx/Nx, Ly/Ny）を自動計算する
3. When Z_grid.typeが"uniform"の場合, the Solver shall Lzから等間隔Z座標を自動生成する
4. When Z_grid.typeが"non-uniform"の場合, the Solver shall テキストファイルからセル界面座標を読み込む
5. When Z_grid.typeが"non-uniform"の場合, the Solver shall ファイルの格子点数がNz+5と一致するか検証し、不一致時はエラー終了する
6. The Solver shall 入力として領域寸法(Lx, Ly)、内部セル数(Nx, Ny, Nz)を受け付ける
7. The Solver shall Z_grid.typeに応じてLzまたはZ座標ファイルパスを受け付ける
8. The Solver shall 3次元計算専用として動作する

#### 用語定義

| 用語 | 定義 |
| ---- | ---- |
| セル | 有限体積法の制御体積。物理量（速度、圧力等）はセル中心に配置 |
| 格子点 | セルの界面（稜線）上の点。隣接セルの境界を定義 |
| 内部セル | 計算対象となる物理領域内のセル |
| ゴーストセル | 境界条件適用のため内部セルの外側に設けた仮想セル |
| セル幅 | 隣接する格子点間の距離（= セルのサイズ） |
| Nx, Ny, Nz | 各軸方向の内部セル数 |

```
格子点（界面）とセルの関係（1次元の例、Nz=3の場合）:

格子点:    ●────●────●────●────●────●────●────●
           1    2    3    4    5    6    7    8
                    ↑              ↑
                  z_min          z_max
                  境界            境界

セル:         [G1] [G2] [ 1] [ 2] [ 3] [G3] [G4]
              ゴースト   内部セル(Nz=3)   ゴースト

- 格子点数: Nz + 5 = 8点
- 総セル数: Nz + 4 = 7セル（内部3 + ゴースト4）
```

#### Z方向格子指定

| Z_grid.type | 必須パラメータ | 動作 |
| ----------- | -------------- | ---- |
| uniform | Lz | Lzから等間隔Z座標を自動生成（Δz = Lz/Nz） |
| non-uniform | Z_grid.file | テキストファイルからセル界面座標を読み込み |

#### Z座標ファイルフォーマット（non-uniform時）
```
<格子点数>
<番号> <座標値>
...
```

- 格子点数はセル界面数（HALO領域含む）を表す。Nz個の物理セルに対して、格子点数 = Nz + 5 となる（全Nz+4セルを区切る界面点数）
- 番号は1オリジン（Julia準拠）。1から格子点数（Nz+5）までの連番
  - 番号1〜2: z_min側HALOセル界面（外挿点）
  - 番号3〜Nz+3: 物理領域セル界面（番号3がz_min境界、番号Nz+3がz_max境界）
  - 番号Nz+4〜Nz+5: z_max側HALOセル界面（外挿点）
- 座標値は`Origin_of_Region[3]`（Z方向基点）からの相対座標[m]で記述する

### Requirement 2: 支配方程式・乱流モデル
**Objective:** As a シミュレーションエンジニア, I want 弱圧縮性Navier-Stokes方程式とLESモデルを使用できる, so that クリーンルーム内の乱流場を解析できる

#### Acceptance Criteria
1. The Solver shall 弱圧縮性を考慮した非圧縮性Navier-Stokes方程式を解く（連続の式はref_eq.pdfの1.27式に準拠）
2. The Solver shall 標準Smagorinskyモデルによる乱流粘性を計算する（詳細はref_eq.pdf pp.13-14に準拠）
3. The Solver shall Smagorinsky定数Csをパラメータとして受け付ける（デフォルト: 0.2）

### Requirement 3: 入力データの無次元化
**Objective:** As a シミュレーションエンジニア, I want 有次元の入力パラメータを無次元化して計算できる, so that 数値的に安定した計算ができる

#### Acceptance Criteria
1. The Solver shall 読み込んだ有次元パラメータを、代表長さ(L₀)・代表速度(U₀)を用いて無次元量に変換する
2. The Solver shall 内部計算を無次元量で実行する
3. The Solver shall SPH出力を常に有次元量で出力する
4. The Solver shall チェックポイントを常に無次元量で保存する
5. The Solver shall history/コンソール出力を無次元量で出力する
6. The Solver shall 可視化画像/断面テキストを有次元量で出力する

#### 無次元化規則

| 物理量 | 無次元化 | 備考 |
| ------ | -------- | ---- |
| 長さ | x* = x / L₀ | L₀: Reference_Length |
| 速度 | u* = u / U₀ | U₀: Reference_Velocity |
| 時間 | t* = t / (L₀/U₀) | 代表時間 = L₀/U₀ |
| 密度 | ρ* = 1 | ρ₀ = 1（一定密度、入力パラメータ不要） |
| 圧力 | p* = p / U₀² | 動圧基準（ρ₀ = 1のため p* = p / U₀²） |
| レイノルズ数 | Re = U₀L₀/ν | ν: Kinematic_Viscosityから計算 |

#### 出力の次元規約

| 出力種別 | 次元 | 備考 |
| -------- | ---- | ---- |
| SPHファイル（vel, prs, vel_avg） | 有次元 | 変換設定は持たない |
| チェックポイント | 無次元 | is_dimensionalは常に0 |
| historyファイル | 無次元 | 常に無次元量で出力 |
| コンソール表示 | 無次元 | 常に無次元量で表示 |
| 可視化画像/断面テキスト | 有次元 | 可視化は有次元で出力 |

### Requirement 4: 空間離散化
**Objective:** As a シミュレーションエンジニア, I want 有限体積法で空間離散化できる, so that 保存則を満たす数値解が得られる

#### Acceptance Criteria
1. The Solver shall 有限体積法・セルセンター変数配置を使用する
2. The Solver shall 対流項にWENO3スキームを適用する
3. The Solver shall 拡散項およびその他の項に2次精度中心差分を適用する
4. The Solver shall Lax-Friedrichsフラックス分割を使用する
5. The Solver shall 非等間隔格子に対応したWENO再構成を実装する

#### WENO3フラックス分割

**Lax–Friedrichs（Rusanov）分割**を使用する：

```
f(u) = f⁺(u) + f⁻(u)
f±(u) = (1/2)(f(u) ± α·u)
```

- α = 各成分の最大速度（x方向: max|u|, y方向: max|v|, z方向: max|w|）
- f⁺：左から右へWENO3再構成
- f⁻：右から左へWENO3再構成（左右反転ステンシル）
- 各成分独立に再構成（Characteristic分解は不要）

#### WENO3境界処理

**必要ゴーストセル数:** 2セル（WENO3は2点ステンシル×2）

| 条件 | ゴーストセル設定 |
| ---- | ---------------- |
| Dirichlet | 境界値を直接代入 |
| Neumann | 内部勾配をコピー |
| 壁（反射） | 法線速度符号反転、他はコピー |

※ フラックス分割後の量に対して境界条件を課す
※ 境界付近では精度は最大2次まで低下（仕様として許容）

#### WENO3非等間隔格子再構成

**前提:**
- セル境界：x_{i±1/2}
- セル幅：Δx_i = x_{i+1/2} - x_{i-1/2}
- セル中心：x_i
- 入力はセル平均値 ū_i
- 再構成点：x_{i+1/2}（左側再構成 u⁻_{i+1/2}）

**候補多項式（線形再構成）:**

セル境界 x_{i+1/2} において：
- Stencil 0：{i-1, i}
- Stencil 1：{i, i+1}

各stencilで局所線形補間：
```
q_k = a_k·ū_{i-1} + b_k·ū_i + c_k·ū_{i+1}
```

**再構成係数（具体式）:**

Stencil 0（{i-1, i}）：
```
a₀ = -Δx_i / (Δx_{i-1} + Δx_i)
b₀ = (Δx_{i-1} + 2·Δx_i) / (Δx_{i-1} + Δx_i)
c₀ = 0
```

Stencil 1（{i, i+1}）：
```
a₁ = 0
b₁ = Δx_{i+1} / (Δx_i + Δx_{i+1})
c₁ = Δx_i / (Δx_i + Δx_{i+1})
```

※ 等間隔（Δx_i = Δx）では q₀ = (3/2)ū_i - (1/2)ū_{i-1}, q₁ = (1/2)ū_i + (1/2)ū_{i+1} となり標準WENO3に一致

**理想重み（3次精度条件）:**

```
d₀ = Δx_i / (Δx_{i-1} + 2·Δx_i + Δx_{i+1})
d₁ = (Δx_{i-1} + Δx_i + Δx_{i+1}) / (Δx_{i-1} + 2·Δx_i + Δx_{i+1})
```
（d₀ + d₁ = 1）

※ 等間隔では d₀ = 1/3, d₁ = 2/3

**滑らかさ指標（非等間隔対応）:**
```
β₀ = ((ū_i - ū_{i-1}) / (x_i - x_{i-1}))² · Δx_ref²
β₁ = ((ū_{i+1} - ū_i) / (x_{i+1} - x_i))² · Δx_ref²
```

**代表セル幅 Δx_ref:**
```
Δx_ref = (1/3)(Δx_{i-1} + Δx_i + Δx_{i+1})
```
- 目的：β_k を無次元化し、格子歪みに対する重みの過敏性を低減

**非線形重み:**
```
α_k = d_k / (ε + β_k)^p
ω_k = α_k / (α₀ + α₁)
```
- p = 2, ε = 10⁻⁶

**最終再構成:**
```
u⁻_{i+1/2} = ω₀·q₀ + ω₁·q₁
```

### Requirement 5: 時間積分
**Objective:** As a シミュレーションエンジニア, I want 時間積分スキームを選択できる, so that 精度と計算コストのバランスを調整できる

#### Acceptance Criteria
1. The Solver shall Euler陽解法による時間積分を実行する（必須）
2. Where RK2が選択された場合, the Solver shall 2次Runge-Kutta法を適用する（オプション）
3. Where RK4が選択された場合, the Solver shall 4次Runge-Kutta法を適用する（オプション）
4. The Solver shall 固定時間刻みを使用する

#### 固定時間刻みの計算

時間刻みはシミュレーション開始時に一度だけ計算し、全ステップで固定値として使用する：

```
Δt* = Co · Δx*_min / U*_ref
```

- Co: クーラン数（入力パラメータ）
- Δx*_min: 最小格子幅（無次元）= min(Δx, Δy, min(Δz))
- U*_ref: 参照速度（無次元）= max(U*_max, 1.0)
  - U*_max = sqrt(u_max² + v_max² + w_max²)
  - 初期化時の安定性のため下限値1.0を設定

### Requirement 6: 圧力-速度分離
**Objective:** As a シミュレーションエンジニア, I want Fractional Step法で圧力と速度を分離できる, so that 非圧縮性条件を満たす速度場が得られる

#### Acceptance Criteria
1. The Solver shall 弱い圧縮性を考慮した連続の式（Eq 1.27）に対応するFractional Step法により圧力-速度分離を行う
2. The Solver shall 対流流束・粘性流束から擬似速度をセルセンターで計算する
3. The Solver shall 擬似速度をセルフェイス（スタガード位置）に内挿する
4. The Solver shall セルフェイス擬似速度からコンパクトステンシルでセルの発散値を計算し、ポアソン方程式のソース項とする
5. The Solver shall 擬似速度場に対して速度境界条件を適用してからフェイス内挿を行う
6. The Solver shall 圧力勾配によりセルフェイス速度を修正する
7. The Solver shall 両側フェイスの圧力勾配の平均によりセルセンター速度を修正する
8. The Solver shall 発散値は常にセルフェイス速度を用いて計算する
9. The Solver shall 壁面境界ではマスクによりNeumann条件を課す（圧力勾配にm=0を掛ける）

#### スタガード配置への内挿

セルセンター変数配置でFractional Step法を適用する場合、速度場と圧力場のデカップリング（チェッカーボード不安定性）を防止するため、以下の手順で発散を計算する：

1. **擬似速度計算**: セルセンターで擬似速度ベクトル (u*, v*, w*) を計算
2. **セルフェイスへの内挿**: 擬似速度を各軸方向のセルフェイスに内挿
   ```
   u*_{i+1/2,j,k} = (u*_{i,j,k} + u*_{i+1,j,k}) / 2
   v*_{i,j+1/2,k} = (v*_{i,j,k} + v*_{i,j+1,k}) / 2
   w*_{i,j,k+1/2} = (w*_{i,j,k} + w*_{i,j,k+1}) / 2
   ```
3. **発散計算（コンパクトステンシル）**: セルフェイス値を用いてセルの発散を計算
   ```
   div_{i,j,k} = (u*_{i+1/2} - u*_{i-1/2})/Δx + (v*_{j+1/2} - v*_{j-1/2})/Δy + (w*_{k+1/2} - w*_{k-1/2})/Δz
   ```
   コンパクトステンシルにより、チェッカーボード不安定性を防止し、保存性を保証する
4. **ポアソン方程式**: この発散値がソース項となる
5. **セルフェイス速度の修正**: 圧力勾配でセルフェイス速度を修正
   ```
   u^{n+1}_{i+1/2} = u*_{i+1/2} - Δt * (p_{i+1} - p_i)/Δx
   ```
6. **セルセンター速度の修正**: 両側フェイスの圧力勾配の平均でセルセンター速度を修正
   ```
   u^{n+1}_i = u*_i - Δt * 0.5 * [(p_i - p_{i-1})/Δx * m_{i-1} + (p_{i+1} - p_i)/Δx * m_{i+1}]
   ```
   ここで m はマスク値（壁面では0、流体では1）で、Neumann条件を課す

※ この内挿処理により、隣接セル間の圧力結合が保証され、チェッカーボードパターンの発生を抑制する

### Requirement 7: 圧力ポアソン方程式ソルバー
**Objective:** As a シミュレーションエンジニア, I want 圧力ポアソン方程式を反復法で解ける, so that 圧力場を効率的に計算できる

#### Acceptance Criteria
1. The Solver shall Red-Black SOR法による反復解法を実装する（必須）
2. Where CG法が選択された場合, the Solver shall 共役勾配法を適用する（オプション）
3. Where BiCGSTAB法が選択された場合, the Solver shall BiCGSTAB法を適用する（オプション）
4. The Solver shall 加速係数、収束判定値、最大反復回数をパラメータとして受け付ける
5. The Solver shall SOR残差を初期残差で正規化して評価する（H2方式）
6. The Solver shall CG/BiCGSTABで前処理付き共役勾配法（Gauss-Seidel 5 sweep）を適用する
7. The Solver shall 収束後に圧力場の空間平均値を計算し、外部境界に圧力基準点（Outflow等）が存在しない場合に限り全セルから平均値を減算する（ドリフト防止）

#### 圧力平均値の引き戻し

圧力境界条件が全面Neumannの場合、ポアソン方程式の解は定数分の不定性を持つ。
時間積分を繰り返すと圧力値がドリフトするため、外部に固定圧力境界（Outflow等）が存在しない場合に限り、各ステップで以下の処理を行う：

```
p_avg = (1/N) Σ p_i      # 流体セルの空間平均
p_i = p_i - p_avg        # 全流体セルから平均値を減算
```

※ 物体セル（mask=0）は平均値計算から除外する
※ Outflow境界など、圧力がDirichletで規定される境界が一つでも存在する場合、この平均値減算は行わない（基準圧力が固定されるため）

#### 残差定義

線形方程式 Ax = b に対して、残差を以下で定義する：

```
SOR残差 = ||r||₂ / Res0
CG/BiCGSTAB残差 = ||r||₂ / ||r0||₂
```

**SOR残差（真の残差方式）**

```
r = b - A p    # 真の残差 (b_val - (ss - dd*pp))
Res0 = ||r||₂（初期反復時）
```

- Res0 = 0 の場合は Res0 = 1 に置き換える

**CG/BiCGSTAB残差**

```
r0 = b - A x0
```

- r0 = 0 の場合は収束済みとして終了

### Requirement 8: 境界条件
**Objective:** As a シミュレーションエンジニア, I want 各種境界条件を設定できる, so that 実際のクリーンルーム環境を模擬できる

#### Acceptance Criteria
1. The Solver shall 外部境界（6面）に対して標準条件（wall, symmetric, periodic, outflow, SlidingWall）および特殊条件（neumann, dirichlet）を速度に適用できる
2. The Solver shall 外部境界の圧力にNeumann条件を適用する
3. The Solver shall 壁面にデフォルトで粘着条件（速度）とNeumann条件（圧力）を適用する
4. When 吹出口・吸込口が指定された時, the Solver shall 座標指定で速度を設定する
5. The Solver shall 領域内の部分境界（矩形/円筒領域）に法線・速度を指定できる
6. The Solver shall 対流流出条件の輸送速度を該当境界面の平均速度とする
7. The Solver shall 周期境界条件を外部境界のペア（x_min/x_max, y_min/y_max, z_min/z_max）に適用できる
8. The Solver shall 対称境界条件で法線方向速度成分をゼロ、接線方向速度成分の勾配をゼロとする

| 種別               | 速度                                            | 圧力    |
| ------------------ | ----------------------------------------------- | ------- |
| 外部境界（標準）    | 壁 (wall) / 対称 (symmetric) / 周期 (periodic) / 対流流出 (outflow) / SlidingWall | Neumann / 周期 / Dirichlet (outflowのみ) |
| 外部境界（特殊）    | Dirichlet / Neumann | Neumann |

※ `outflow` 境界では、圧力の安定性を確保するため基準圧力として `prs = 0.0` (Dirichlet) を適用する。

※ `wall`, `symmetric`, `SlidingWall` の場合、ゴーストセルのマスク値は `0`（固体）に設定される。これに伴い、拡散項の計算ではマスク値に基づき境界での勾配がデフォルトでゼロ（Neumann）となる。
※ `wall` および `SlidingWall` においては、物理的なノンスリップ条件を満足するため、拡散項の計算後に壁面せん断流束（例：Z+面の場合 $2 \nu_{eff} (U_w - u) / \Delta z^2$）を明示的に加算して修正を行う。
※ `SlidingWall` は境界速度を指定する壁面（例: キャビティ流れの上壁）に使用する。指定値は `value` で与える。
※ JSON で指定されるすべての文字列パラメータ（境界条件名、時間スキーム、ソルバー名等）は、**読み取り時に大文字小文字を区別せず評価される**（Case-Insensitive）。

※ Neumann条件は境界面の法線方向勾配（∂φ/∂n）を指定する

#### 周期境界条件

対向する境界面（例: y_min/y_max）間で周期性を持たせる：

```
Ghost(min側) ← Inner(max側)
Ghost(max側) ← Inner(min側)
```

- 対向面の両方がPeriodicである必要がある（片方のみはエラー）
- 速度・圧力の両方にゴーストセル交換を適用
- 2Dシミュレーション（XZ平面）では y_min/y_max に適用

#### 対流流出条件

各速度成分φに対して以下の式を適用する：

```
∂φ/∂t + Uc · ∂φ/∂n = 0
```

境界面の法線方向に応じて輸送速度Ucと離散化が異なる：

| 境界面 | 法線方向 | 輸送速度Uc | 離散化（1次風上） |
| ------ | -------- | ---------- | ----------------- |
| x_min / x_max | x | 境界面の平均u | φ_new = φ - Uc·Δt·(φ - φ_neighbor)/Δx |
| y_min / y_max | y | 境界面の平均v | φ_new = φ - Uc·Δt·(φ - φ_neighbor)/Δy |
| z_min / z_max | z | 境界面の平均w | φ_new = φ - Uc·Δt·(φ - φ_neighbor)/Δz |

- φ: 速度成分（u, v, w のいずれか）
- φ_neighbor: 境界に隣接する内部セルの値
- Δx, Δy, Δz: セル幅（非等間隔格子では局所セル幅を使用）

※ 壁面（粘着条件）は `"velocity": "dirichlet", "value": [0.0, 0.0, 0.0]` で指定


#### 対称境界条件

対称境界条件は、境界面に対して法線方向の速度成分をゼロ、接線方向の速度成分の勾配をゼロとする条件です。

**速度境界条件:**

境界面の法線方向をn、接線方向をt1, t2とすると：

\`\`\`
u_n = 0           (法線方向速度成分 = 0)
∂u_t1/∂n = 0      (接線方向速度成分の法線勾配 = 0)
∂u_t2/∂n = 0      (接線方向速度成分の法線勾配 = 0)
\`\`\`

**各境界面での実装:**

| 境界面 | 法線方向 | 条件 |
|--------|----------|------|
| x_min/x_max | x | u=0, ∂v/∂x=0, ∂w/∂x=0 |
| y_min/y_max | y | v=0, ∂u/∂y=0, ∂w/∂y=0 |
| z_min/z_max | z | w=0, ∂u/∂z=0, ∂v/∂z=0 |

**ゴーストセルへの適用:**

法線方向速度成分は符号反転、接線方向速度成分はミラーリング：

\`\`\`julia
# x_min境界の例（i=3が境界）
u[2, j, k] = -u[3, j, k]  # 法線方向: 符号反転
u[1, j, k] = -u[4, j, k]
v[2, j, k] = v[3, j, k]   # 接線方向: ミラーリング
v[1, j, k] = v[4, j, k]
w[2, j, k] = w[3, j, k]
w[1, j, k] = w[4, j, k]
\`\`\`

**圧力境界条件:**

対称境界では圧力の法線勾配をゼロとする（Neumann条件）：
\`\`\`
∂p/∂n = 0
\`\`\`

**用途:**

- 対称性を持つ流れ場の計算領域を半分にすることで計算コストを削減
- 例: チャネル流れの中心面、対称な物体周りの流れ

#### 境界条件の優先順位

セルが複数の境界条件の対象となる場合、以下の優先順位で適用する：

| 優先度 | 種別 | 説明 |
| ------ | ---- | ---- |
| 1（最高） | 内部境界 | 吹出口・吸込口、領域内の部分境界 |
| 2 | 物体 | ボクセル法で表現された障害物 |
| 3（最低） | 外部境界 | 計算領域の6面 |

### Requirement 9: 物体表現
**Objective:** As a シミュレーションエンジニア, I want 計算領域内に障害物を配置できる, so that クリーンルーム内の設備を模擬できる

#### Acceptance Criteria
1. The Solver shall ボクセル法（マスク法）により物体を表現する
2. The Solver shall セルマスク値（流体=1、物体=0）によりfluxを制御する
3. The Solver shall Geometry JSONファイルから物体形状定義を読み込む
4. The Solver shall 直方体（box）、円筒（cylinder）、球（sphere）の基本形状をサポートする
5. When 基本形状が定義された時, the Solver shall 形状内部に含まれるセルを物体セル（マスク値=0）として設定する
6. The Solver shall 任意の数の物体を配置できる
7. The Solver shall 物体内部の速度を指定速度Vsに設定する（デフォルト: [0.0, 0.0, 0.0]）
8. The Solver shall 物体内部の圧力をポアソン方程式により計算する（セルフェイス値から自然に決定）

#### Geometry JSON仕様

```json
{
  "objects": [
    {
      "name": "equipment_1",
      "type": "box",                       // (box | cylinder | sphere)
      "min": [1.0, 1.0, 0.0],              // [m] 最小座標
      "max": [2.0, 2.0, 1.5],              // [m] 最大座標
      "velocity": [0.0, 0.0, 0.0]          // [m/s] 物体内部速度Vs（デフォルト: [0.0, 0.0, 0.0]）
    },
    {
      "name": "pillar_1",
      "type": "cylinder",
      "center": [3.0, 3.0, 0.0],           // [m] 底面中心座標
      "radius": 0.3,                       // [m] 半径
      "height": 2.5,                       // [m] 高さ
      "axis": "z",                         // (x | y | z) 軸方向
      "velocity": [0.0, 0.0, 0.0]          // [m/s] 物体内部速度Vs（デフォルト: [0.0, 0.0, 0.0]）
    },
    {
      "name": "obstacle_1",
      "type": "sphere",
      "center": [2.5, 2.5, 1.5],           // [m] 中心座標
      "radius": 0.5,                       // [m] 半径
      "velocity": [0.0, 0.0, 0.0]          // [m/s] 物体内部速度Vs（デフォルト: [0.0, 0.0, 0.0]）
    }
  ]
}
```

#### 物体内部の処理

| 物理量 | 処理 | 備考 |
| ------ | ---- | ---- |
| 速度 | 指定速度Vsを設定 | JSONの`velocity`で指定、省略時は[0.0, 0.0, 0.0] |
| 圧力 | ポアソン方程式で計算 | セルフェイス値から右辺項が決まり自然に決定 |

※ 物体内部の速度指定は完全な不動壁を前提とする（移動物体・回転物体は非対応）

### Requirement 10: HALO領域
**Objective:** As a シミュレーションエンジニア, I want HALO領域を確保できる, so that 将来の分散並列計算に対応できる

#### Acceptance Criteria
1. The Solver shall Z方向の外部境界（z_min側およびz_max側の両方）にステンシル幅分のHALO領域を確保する

※ XY方向は等間隔格子であり、HALO領域は不要（Z方向のみ非等間隔格子をサポート）

### Requirement 11: 入出力
**Objective:** As a シミュレーションエンジニア, I want 各種ファイル形式で入出力できる, so that パラメータ設定と可視化ツールとの連携ができる

#### Acceptance Criteria
1. The Solver shall 計算パラメータをJSON形式で読み込む
2. The Solver shall 境界条件をJSON形式で読み込む
3. When Z_grid.typeが"non-uniform"の場合, the Solver shall Z方向格子座標をテキスト形式で読み込む
4. When start="initial"の場合, the Solver shall 初期条件（一様速度ベクトル）をJSON形式で読み込む
5. The Solver shall 速度場および圧力場をV-Isio SPH形式（バイナリ、単精度、Fortran unformatted互換）で有次元量として出力する
6. The Solver shall チェックポイントデータを無次元のヘッダ付き倍精度バイナリ形式で指定時間間隔毎に保存する（ヘッダ: 格子数、ステップ数、時刻、次元フラグ=0、代表量）
7. When start="restart"の場合, the Solver shall Restart.fileで指定されたチェックポイントファイルから計算を再開する
8. When start="restart"かつチェックポイントの次元フラグが0以外の場合, the Solver shall エラーメッセージを出力し計算を停止する
9. When 計算パラメータが読み込まれた時, the Solver shall CFL条件（クーラン数）と拡散数による安定条件をチェックする
10. If 安定条件を満たさない場合, then the Solver shall 警告を表示し計算を停止する（dry_run時は警告のみ）
11. The Solver shall 指定間隔で計算モニター情報をコンソールに表示する
12. The Solver shall historyファイルに各ステップのモニター値をテキスト形式で出力する
13. The Solver shall 速度場の時間平均値を指定間隔で可視化フォーマット（SPH形式）で出力する
14. The Solver shall 時間平均値をインクリメンタル方式で計算する
15. If 入力ファイルが不正な場合, then the Solver shall エラーメッセージを出力し計算を停止する
16. If 計算中に発散が検出された場合, then the Solver shall エラーメッセージを出力し計算を停止する
17. When dry_runが有効な場合, the Solver shall 必要メモリ量（作業配列の合計）を計算し表示する
18. When dry_runが有効な場合, the Solver shall システムの利用可能メモリと比較し実行可否を判定する

#### エラーハンドリング

| エラー種別 | 検出条件 | 動作 |
| ---------- | -------- | ---- |
| 入力ファイル不正 | JSON構文エラー、必須パラメータ欠落、型不一致 | エラーメッセージ出力、停止 |
| 安定条件違反 | 拡散数 D ≥ 0.5 | 警告出力、停止（dry_run時は警告のみ） |
| メモリ不足 | 必要メモリ > 利用可能メモリ | 警告出力（dry_run時のみチェック） |
| 発散検出 | divMax > divMax_threshold、NaN/Inf検出 | エラーメッセージ出力、停止 |
| ファイルI/Oエラー | 読込/書込失敗 | エラーメッセージ出力、停止 |

#### 計算モニター項目

コンソール表示およびhistoryファイル出力に共通:

| 項目 | 説明 |
| ---- | ---- |
| step | タイムステップ数（左詰め、可変幅パディング） |
| time | 計算時刻（%.6e） |
| Umax | 速度の最大値（%.4e） |
| divMax | 速度の発散値の最大値（%.4e） |
| dU | 速度変動量（%.4e）定常状態判定用 |
| ItrP | 圧力反復回数（整数） |
| ResP | 圧力反復残差（%.5e） |

※ divMax = max|∂u/∂x + ∂v/∂y + ∂w/∂z|（全セルの絶対値最大）
※ dU = sqrt(Σ((u^n+1 - u^n)² + (v^n+1 - v^n)² + (w^n+1 - w^n)²))（全内部セルの速度変動のL2ノルム）

#### historyファイル仕様

- 形式: テキスト（スペース区切り、固定幅カラム）
- ファイル名: `history.txt`
- 出力場所: `dirname(param_file)/output/history.txt`
- 出力間隔: `Intervals.history` で指定（ステップ単位）
- フォーマット:
  - Header: 右詰め（stepのみ左詰め）
  - Data: 右詰め（stepのみ左詰め）
  - 位置合わせされたカラム出力

#### condition.txtファイル仕様

シミュレーション開始前に計算条件を出力する：

- ファイル名: `condition.txt`
- 出力場所: `dirname(param_file)/output/condition.txt`
- 出力項目:
  - 物理パラメータ（L0, U0, ν, Re, T0）
  - 格子パラメータ（格子数、ドメインサイズ、セルサイズ）
  - 時間積分（スキーム、Co、Δt*、Δt、最大ステップ、総計算時間）
  - ポアソンソルバー（ソルバー種別、収束判定値、最大反復数）
  - 出力間隔

#### 速度の時間平均値

- 対象: 速度場（u, v, w）のみ（圧力は平均化対象外）
- 平均化開始時刻: `Start_time_for_averaging` で指定（JSON入力は有次元[sec]、読み込み後は無次元化して保持: t* = t / (L₀/U₀)）
- 出力間隔: `Intervals.averaged_file` で指定（ステップ単位）
- 可視化用出力: SPH形式（単精度）
- 計算方式: インクリメンタル平均（Welfordの方法）

```
# インクリメンタル平均の更新式
n = n + 1
mean_new = mean_old + (x - mean_old) / n
```

※ メモリ効率が良く、数値的に安定

#### SPHファイル仕様（単精度・ベクトル）

出力は常に有次元量で行う。

エンディアン: **リトルエンディアン**

| レコード   | 内容                         | 備考                                    |
| ---------- | ---------------------------- | --------------------------------------- |
| データ属性 | svType=2, dType=1            |                                         |
| サイズ     | IMAX, JMAX, KMAX (Int32)     |                                         |
| 原点座標   | XORG, YORG, ZORG (Float32)   | **最初のセルの左端（フェイス位置）**    |
| ピッチ     | DX, DY, DZ (Float32)         | 不等間隔でも等間隔値(Lx/Nx等)を格納     |
| 時刻       | step (Int32), time (Float32) |                                         |
| データ     | 速度成分 (Float32)           | u, v, w の順に配列を記述               |

**座標系の注意:**
- 原点座標（XORG, YORG, ZORG）は最初のセルの**左端（フェイス位置）**を示す
- セルセンター座標は `x_center[i] = XORG + (i - 0.5) * DX` で計算
- これにより、最後のセルセンターは `XORG + (N - 0.5) * DX` となり、計算領域端（壁面位置）には達しない

#### SPHファイル仕様（単精度・スカラー）

出力は常に有次元量で行う。

エンディアン: **リトルエンディアン**

| レコード   | 内容                         | 備考                                    |
| ---------- | ---------------------------- | --------------------------------------- |
| データ属性 | svType=1, dType=1            | svType=1: スカラー                      |
| サイズ     | IMAX, JMAX, KMAX (Int32)     |                                         |
| 原点座標   | XORG, YORG, ZORG (Float32)   | **最初のセルの左端（フェイス位置）**    |
| ピッチ     | DX, DY, DZ (Float32)         | 不等間隔でも等間隔値(Lx/Nx等)を格納     |
| 時刻       | step (Int32), time (Float32) |                                         |
| データ     | 圧力 (Float32)               | スカラー配列を1つ記述                   |

※ 圧力SPHファイル（prs_#######.sph）はこの形式で出力する

#### ファイル命名規則

| ファイル種別 | 命名規則 | 例 | 精度 |
| ------------ | -------- | -- | ---- |
| 瞬時値速度 | `vel_#######.sph` | vel_0000100.sph | 単精度 |
| 瞬時値圧力 | `prs_#######.sph` | prs_0000100.sph | 単精度 |
| 平均値速度（可視化用） | `vel_avg_#######.sph` | vel_avg_0000100.sph | 単精度 |
| チェックポイント | `checkpoint_#######.bin` | checkpoint_0010000.bin | 倍精度 |

※ `#######` はタイムステップ数（7桁、ゼロ埋め）

#### チェックポイントファイル仕様（倍精度バイナリ）

エンディアン: **リトルエンディアン**

| 項目 | 型 | 説明 |
| ---- | -- | ---- |
| **ヘッダ部** | | |
| Nx, Ny, Nz | Int32 × 3 | 内部セル数 |
| step | Int32 | タイムステップ数 |
| time | Float64 | 計算時刻 |
| is_dimensional | Int32 | 0固定（無次元のみ） |
| Reference_Length | Float64 | 代表長さ L₀ [m] |
| Reference_Velocity | Float64 | 代表速度 U₀ [m/s] |
| **データ部** | | |
| u | Float64[Nx+4, Ny+4, Nz+4] | 速度u成分（ゴーストセル含む） |
| v | Float64[Nx+4, Ny+4, Nz+4] | 速度v成分（ゴーストセル含む） |
| w | Float64[Nx+4, Ny+4, Nz+4] | 速度w成分（ゴーストセル含む） |
| p | Float64[Nx+4, Ny+4, Nz+4] | 圧力（ゴーストセル含む） |

※ 配列サイズ(Nx+4, Ny+4, Nz+4)はWENO3の5点ステンシルに対応するゴーストセル（各軸両側2セル）を含む

#### 次元フラグの扱い

- is_dimensionalは常に0（無次元）
- 0以外の場合はフォーマット不正としてエラー終了する

※ リスタート時のJSON設定（Reference_Length, Reference_Velocity）とチェックポイントヘッダの値が一致することを確認し、不一致の場合は警告を出力する

#### 初期条件とリスタート

| start | 動作 | 必須パラメータ |
| ----- | ---- | -------------- |
| initial | 一様初期条件で新規開始 | Initial_Condition.velocity, Initial_Condition.pressure |
| restart | チェックポイントから再開 | Restart.file |

※ リスタート時はチェックポイントファイルに保存されたステップ数・時刻から計算を継続

#### 計算パラメータJSON仕様

```json
{
  "dry_run": "no",                        // (yes | no) ドライラン
  "start": "initial",                     // (initial | restart) 開始モード
  "Max_step": 100000,                     // 計算ステップ数
  "Reference_Length": 0.5,                // [m] 代表長さ
  "Reference_Velocity": 0.3,              // [m/s] 代表速度
  "Kinematic_Viscosity": 1.5e-5,          // [m²/s] 動粘性係数
  "Smagorinsky_Constant": 0.2,            // [-] Smagorinsky定数Cs（デフォルト: 0.2）
  "Origin_of_Region": [-2.5, -2.5, 0.0],  // [m] 計算領域基点（ガイドセルを含まない物理領域の各軸最小値）
  "Domain": {
    "Lx": 5.0,                            // [m] X方向領域長さ
    "Ly": 5.0,                            // [m] Y方向領域長さ
    "Nx": 100,                            // [-] X方向内部セル数
    "Ny": 100,                            // [-] Y方向内部セル数
    "Nz": 50                              // [-] Z方向内部セル数
  },
  "Z_grid": {
    "type": "uniform",                    // (uniform | non-uniform) Z方向格子タイプ
    "Lz": 2.5,                            // [m] Z方向領域長さ（uniform時に有効）
    "file": "z_grid.txt"                  // Z座標ファイル（non-uniform時に有効）
  },
  "Courant_number": 0.2,                  // [-] クーラン数（Δt自動計算用）
  "Intervals": {
    "display": 100,                       // コンソール表示間隔[step]
    "history": 1,                         // 履歴出力間隔[step]
    "Instantaneous_file": 4000,           // 瞬時値ファイル出力間隔[step]
    "averaged_file": 4000,                // 平均値ファイル出力間隔[step]
    "checkpoint": 10000                   // チェックポイント出力間隔[step]
  },
  "Visualization": {
    "interval": 1000,                     // 可視化出力間隔[step]
    "plane": "xy",                        // (xy | xz | yz) 可視化断面
    "plane_index": 15,                    // 断面インデックス（指定軸に垂直な断面位置、1オリジン）
    "variables": ["velocity", "pressure"], // 可視化対象変数
    "output_format": "png",               // (png | svg) 出力画像形式
    "output_dir": "viz",                  // 出力ディレクトリ（output配下）
    "vector_enabled": false,              // ベクトル矢印表示
    "vector_skip": 1,                     // ベクトル間引き数
    "text_output": false                  // 断面テキスト出力
  },
  "Poisson_parameter": {
    "solver": "RedBlackSOR",              // (RedBlackSOR | CG | BiCGSTAB) ソルバー種別
    "coef_acceleration": 0.9,             // RedBlackSOR加速係数（RedBlackSOR時のみ有効）
    "convergence_criteria": 1.0e-3,       // 収束判定値
    "Iteration_max": 100,                 // 最大反復回数
    "on_divergence": "WarnContinue"       // (WarnContinue | Abort) 収束失敗時の動作
  },
  "Start_time_for_averaging": 0.0,        // [sec] 平均化開始時刻
  "Time_Integration_Scheme": "Euler",     // (Euler | RK2 | RK4) 時間積分スキーム
  "divMax_threshold": 1.0e-3,             // [-] 発散検出閾値
  "Initial_Condition": {
    "velocity": [0.0, 0.0, 0.0],          // [m/s] 初期速度ベクトル（一様）
    "pressure": 0.0                        // [Pa] 初期圧力値（一様）
  },
  "Restart": {
    "file": "checkpoint_0010000.bin"      // リスタートファイル（start="restart"時に有効）
  }
}
```

#### 時間刻みと安定条件

**時間刻みの決定:**
- 無次元時間刻みΔt*はクーラン数（Courant_number）から自動計算する
- Δt* = Co · Δx*_min / U*_ref （Δx*_min: 最小セル幅、U*_ref: 参照速度）
- U*_ref = max(U*_max, 1.0) （U*_max: 場の最大速度、無次元代表速度=1.0）
- 3次元では各方向の最小セル幅で判定

※ 初期静止場（U*_max=0）や低速場では無次元代表速度（=1.0）を下限として使用し、ゼロ除算を回避する

**安定条件チェック:**
- 拡散数 D = ν*·Δt*/Δx*² を計算し、D < 0.5 を満たすことを確認
- If 拡散数条件を満たさない場合, then 警告を表示し計算を停止する（dry_run時は警告のみ）

| 条件 | 式 | 用途 |
| ---- | -- | ---- |
| CFL条件（対流） | Co = U·Δt/Δx | Δt決定 |
| 拡散数条件（粘性） | D = ν·Δt/Δx² < 0.5 | チェックのみ |

#### ドライラン時のメモリチェック

dry_run有効時に以下を実行する：

1. **必要メモリ量の計算**
   - 全作業配列のサイズを合計（Float64: 8バイト/要素）
   - 配列サイズ: (Nx+4) × (Ny+4) × (Nz+4)（WENO3の5点ステンシル対応、各軸両側2セルのゴースト）
   - 主要配列: u, v, w, p, mask, および各種作業配列（計15配列）

2. **システム情報の取得**
   - 利用可能メモリ（Sys.free_memory()）
   - 総メモリ（Sys.total_memory()）

3. **実行可否判定**
   - 必要メモリ < 利用可能メモリ: 実行可能
   - 必要メモリ ≥ 利用可能メモリ: 警告出力

4. **表示項目**
   - 必要メモリ量 [GB]
   - 利用可能メモリ [GB]
   - 総メモリ [GB]
   - 実行可否

#### 境界条件JSON仕様

```json
{
  "external_boundaries": {
    "x_min": { "velocity": "dirichlet", "value": [0.3, 0.0, 0.0] },
    "x_max": { "velocity": "outflow" },
    "y_min": { "velocity": "periodic" },
    "y_max": { "velocity": "periodic" },
    "z_min": { "velocity": "dirichlet", "value": [0.0, 0.0, 0.0] },
    "z_max": { "velocity": "dirichlet", "value": [0.0, 0.0, 0.0] }
  },
  "inlets": [
    {
      "type": "rectangular",
      "position": [0.0, 0.0, 2.5],         // [m] 中心座標
      "size": [0.5, 0.5],                  // [m] XY方向サイズ
      "normal": [0, 0, -1],                // 法線方向
      "velocity": [0.0, 0.0, -0.5]         // [m/s] 速度ベクトル
    }
  ],
  "outlets": [
    {
      "type": "rectangular",
      "position": [0.0, 0.0, 0.0],
      "size": [0.3, 0.3],
      "normal": [0, 0, -1],
      "condition": "outflow",              // (outflow | dirichlet) 境界条件種別
      "velocity": [0.0, 0.0, -0.3]         // [m/s] 速度ベクトル（dirichlet時に必須）
    }
  ],
  "internal_boundaries": [
    {
      "type": "rectangular",               // (rectangular | cylindrical)
      "region": {
        "min": [1.0, 1.0, 0.5],            // [m] 領域最小座標
        "max": [1.5, 1.5, 1.0]             // [m] 領域最大座標
      },
      "normal": [1, 0, 0],                 // 法線方向
      "velocity": [0.1, 0.0, 0.0]          // [m/s] 速度ベクトル
    },
    {
      "type": "cylindrical",
      "center": [2.0, 2.0, 1.0],           // [m] 中心座標
      "radius": 0.2,                       // [m] 半径
      "height": 0.5,                       // [m] 高さ
      "axis": "z",                         // 軸方向
      "normal": [0, 0, 1],
      "velocity": [0.0, 0.0, 0.2]
    }
  ]
}
```

#### 吸込口（outlets）の境界条件

| condition | 動作 | velocity |
| --------- | ---- | -------- |
| dirichlet | 指定速度を適用 | 必須（速度ベクトルを指定） |
| outflow | 対流流出条件を適用 | 不要（該当断面の平均速度を輸送速度として使用） |

※ `outlets`は`external_boundaries`とは独立した指定体系であり、同一セルに対して両方が指定された場合は`outlets`が優先される（境界条件の優先順位を参照）

### Requirement 12: 可視化機能
**Objective:** As a シミュレーションエンジニア, I want 計算結果を可視化できる, so that 流れ場を確認・解析できる

#### Acceptance Criteria
1. The Solver shall Julia内で可視化機能を実装する（CairoMakieを使用）
2. The Solver shall Package Extensions機能を用いてCairoMakieへの依存を管理し、非可視化実行時のロード時間を最小化する
2. The Solver shall 可視化処理を倍精度データのまま実行する
3. The Solver shall 任意断面（XY, YZ, XZ平面）を抽出・表示する
4. The Solver shall 速度場・圧力場のコンター表示を行う
5. Where ベクトル表示が有効化された場合, the Solver shall ベクトル表示を行う（オプション）
6. The Solver shall 可視化結果を画像ファイル（PNG等）として出力する
7. Where テキスト出力が指定された場合, the Solver shall 断面データをテキストファイル（数値）として出力する

#### 可視化パラメータ（Visualizationセクション）

| パラメータ | 型 | 説明 |
| ---------- | -- | ---- |
| interval | Int | 可視化出力間隔[step] |
| plane | String | 可視化断面（xy, xz, yz） |
| plane_index | Int | 断面インデックス（指定軸に垂直な断面位置、1オリジン） |
| variables | Array[String] | 可視化対象変数（velocity, pressure） |
| output_format | String | 出力画像形式（png, svg） |
| output_dir | String | 出力ディレクトリ（output配下） |
| vector_enabled | Bool | ベクトル矢印表示の有効化 |
| vector_skip | Int | ベクトル間引き数 |
| text_output | Bool | 断面テキスト出力の有効化 |

#### 可視化出力形式

| 出力形式 | 内容 | 備考 |
| -------- | ---- | ---- |
| 画像ファイル | PNG等の画像形式 | デフォルト出力 |
| テキストファイル | 断面の数値データ | オプション、座標と物理量をスペース区切り、有次元で出力 |



---

## Non-Functional Requirements

### NFR-1: 開発言語
The Solver shall Juliaで実装される

### NFR-2: 拡張性
The Solver shall 将来の分散並列計算への拡張性を考慮したHALO領域設計を持つ

### NFR-4: 並列化方針
The Solver shall 最初は逐次計算で開発し、後にスレッド並列化を実施する

### NFR-3: モジュール性
The Solver shall スキーム・ソルバーをモジュール化し差替え可能とする

### NFR-5: 環境管理
The Solver shall Project.tomlとManifest.tomlによりJulia実行環境を管理する

---

## Out of Scope（参考情報）

| 項目       | 内容                                     |
| ---------- | ---------------------------------------- |
| 目的       | クリーンルーム内微小粒子の発生源特定     |
| 後段解析   | Backward FTLE + Backward LPT             |
| 検証問題   | 3Dキャビティフロー、バックステップ流（Re数・格子数はJSONパラメータで指定） |
| 実装優先度 | 基本スキーム先行、高度なものはオプション |

---

## Decisions Summary

| 項目               | 決定内容                           |
| ------------------ | ---------------------------------- |
| 開発言語           | Julia                              |
| 次元               | 3D専用                             |
| Z格子ファイル      | テキスト（格子数+番号/座標ペア）   |
| SPH出力            | バイナリ、単精度、ピッチは等間隔値 |
| チェックポイント   | バイナリ、倍精度（リスタート用）   |
| バイナリエンディアン | リトルエンディアン                 |
| 可視化精度         | 倍精度（内部データをそのまま使用） |
| 可視化ライブラリ   | 後日選定                           |

---

## Implementation Reference

参考実装: `/Users/Daily/Development/H2/src/`

### 格子・配列構造
- 境界セル付き配列: `(MX, MY, MZ) = (NX+2, NY+2, NZ+2)`
- ゴーストセル: インデックス1と最大値、物理領域は2〜size-1
- Z方向非一様格子: 面座標(Z)、セル中心(ZC)、セル幅(ΔZ)を別配列で管理

### データ構造パターン
```julia
struct WorkBuffers
    u::Array{Float64,3}       # 速度u成分
    v::Array{Float64,3}       # 速度v成分
    w::Array{Float64,3}       # 速度w成分
    p::Array{Float64,3}       # 圧力
    mask::Array{Float64,3}    # 物体マスク (1.0=流体, 0.0=物体)
    # ... ソルバー作業配列
end
```

### 幾何形状の格子転写
1. **IDマップ方式**: `ID::Array{UInt8,3}`で物体IDを管理
2. **判定関数**:
   - `is_included_rect`: 直方体（体積50%以上オーバーラップ）
   - `is_included_cyl`: 円筒（サンプリング点法）
   - `is_included_sph`: 球（頂点判定+サンプリング）
3. **バウンディングボックス最適化**: 格子探索範囲を限定

### 並列化パターン（将来実装）
- 初期開発: 逐次計算
- 将来: `FLoops.jl`によるスレッド並列
- バックエンド抽象化: sequential/thread切替
- Red-Black SORでの`@simd`活用
