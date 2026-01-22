# Implementation Plan

## Tasks

- [x] 1. プロジェクト基盤構築
- [x] 1.1 Juliaパッケージ構成とモジュール骨格の作成
  - Project.toml/Manifest.tomlでFLoops.jl, JSON3.jl, LinearAlgebra.jl依存を定義
  - メインモジュール(CleanroomSolver)と各サブモジュールのファイル骨格を作成
  - CairoMakieをweakdepsとして設定し、Package Extensions構成を準備
  - _Requirements: NFR-1, NFR-5_

- [x] 1.2 (P) 共通定数・型・ユーティリティの実装
  - 無次元化パラメータ構造体(DimensionParams)の定義
  - 並列バックエンド切り替え関数(get_backend)の実装
  - Float64精度での数値計算用ヘルパー関数
  - テスト: 無次元化/有次元化の往復変換テスト
  - _Requirements: 3.1, 3.2_

- [x] 2. 格子生成機能
- [x] 2.1 格子構成データ構造の実装
  - GridConfig構造体（Nx, Ny, Nz, Lx, Ly, origin, z_type, z_file）の定義
  - GridData構造体（mx, my, mz, dx, dy, dz, x, y, z_face, z_center）の定義
  - 配列サイズ計算（WENO3の5点ステンシル対応、各軸両側2セルのゴースト）
  - _Requirements: 1.1, 1.6, 1.8, 10.1_

- [x] 2.2 XY方向等間隔格子生成の実装
  - XY方向の等間隔格子座標を自動計算（Lx/Nx, Ly/Ny）
  - セル中心座標配列の生成
  - ゴーストセル領域を含む配列サイズ管理
  - _Requirements: 1.2, 1.6_

- [x] 2.3 Z方向格子生成の実装
  - uniform時：Lzから等間隔Z座標を自動生成（Δz = Lz/Nz）
  - non-uniform時：テキストファイルからセル界面座標を読み込み
  - 格子点数検証（Nz+5点と一致するか確認、不一致時エラー）
  - Z座標のセル中心・セル幅配列の計算
  - _Requirements: 1.3, 1.4, 1.5, 1.7_

- [x] 2.4 格子の無次元化処理
  - 代表長さL₀による格子座標の無次元化
  - generate_grid関数の完成（GridConfig + DimensionParams → GridData）
  - _Requirements: 3.1_

- [ ] 2.5 格子生成の検証（格子断面可視化）
  - XY/XZ/YZ断面での格子線描画による目視検証
  - 等間隔格子（XY方向）のセル幅均一性確認
  - 非等間隔格子（Z方向）のセル界面座標・セル幅分布の確認
  - ゴーストセル領域の正しい配置確認
  - _Requirements: 1.1, 1.2, 1.3, 1.4_

- [ ] 3. 物理量配列管理
- [x] 3.1 CFDBuffers構造体の実装
  - 速度場配列（u, v, w）の定義と初期化
  - 圧力場配列（p）の定義と初期化
  - 擬似速度配列（u_star, v_star, w_star）の定義
  - マスク配列（mask: 流体=1, 物体=0）の定義
  - _Requirements: 3.2, 9.2_

- [x] 3.2 時間平均・作業配列の実装
  - 時間平均配列（u_avg, v_avg, w_avg, avg_count）の定義
  - ポアソンソルバー作業配列（rhs）の定義
  - フラックス作業配列（flux_u, flux_v, flux_w）の定義
  - 乱流粘性配列（nu_t, nu_eff）の定義
  - セルフェイス速度配列（u_face_x, v_face_y, w_face_z）の定義
  - _Requirements: 11.13, 11.14_

- [x] 3.3 (P) オプション作業配列の実装
  - KrylovBuffers構造体（CG用: r, p, q / BiCGSTAB用: r0, s, t追加）の定義
  - RKBuffers構造体（RK2: u_rk1のみ / RK4: u_rk1〜u_rk4）の定義
  - ソルバー/スキーム選択に応じた動的確保
  - _Requirements: 5.2, 5.3, 7.2, 7.3_

- [x] 3.4 時間平均更新機能の実装
  - Welford法によるインクリメンタル平均計算
  - 平均化開始時刻判定（Start_time_for_averaging）
  - update_time_average!関数の実装
  - テスト: 既知の数列に対するWelford法平均値の精度検証
  - _Requirements: 11.13, 11.14_

- [x] 3.5 メモリ見積もり機能の実装
  - estimate_memory_size関数（時間スキーム・ソルバー依存の配列数計算）
  - check_memory_availability関数（Sys.free_memory()との比較）
  - dry_run時の必要メモリ量表示
  - _Requirements: 11.17, 11.18_

- [x] 4. 入力読込機能
- [x] 4.1 計算パラメータJSON読込の実装
  - SimulationParams構造体の定義（全パラメータを包含）
  - load_parameters関数（JSON解析、必須パラメータ検証）
  - 有次元パラメータの無次元化処理
  - _Requirements: 11.1, 3.1_

- [x] 4.2 境界条件JSON読込の実装
  - ExternalBC, InletOutlet, InternalBoundary構造体の定義
  - BoundaryConditionSet構造体の定義
  - load_boundary_conditions関数の実装
  - _Requirements: 11.2, 8.1, 8.4, 8.5_

- [x] 4.3 初期条件・リスタート設定の実装
  - InitialCondition構造体（velocity, pressure）の定義
  - IntervalConfig構造体（各種出力間隔）の定義
  - リスタートファイルパス設定の読込
  - _Requirements: 11.4, 11.7_

- [x] 4.4 安定条件チェックの実装
  - CFL条件に基づく時間刻み計算（Δt* = Co · Δx*_min / U*_ref）
  - 拡散数条件チェック（D = ν*·Δt*/Δx*² < 0.5）
  - 違反時の警告出力と停止処理（dry_run時は警告のみ）
  - _Requirements: 11.9, 11.10_

- [ ] 4.5 入力読込機能のテスト
  - テスト用JSON設定ファイルの作成（計算パラメータ、境界条件）
  - 正常系: 全パラメータ指定時の読込テスト
  - 異常系: 必須パラメータ欠落時のエラー検出テスト
  - 異常系: 安定条件違反時の警告出力テスト
  - _Requirements: 11.1, 11.2, 11.9, 11.10, 11.15_

- [x] 5. 物体形状・マスク生成
- [x] 5.1 Geometry JSON読込の実装
  - GeometryObject構造体（name, type, params, velocity）の定義
  - load_geometry関数（JSON読込、形状リスト返却）
  - 直方体/円筒/球の基本形状パラメータ解析
  - _Requirements: 9.3, 9.4_

- [x] 5.2 形状判定関数の実装
  - 直方体（box）の内外判定
  - 円筒（cylinder）の内外判定（軸方向対応）
  - 球（sphere）の内外判定
  - セル中心座標に基づく判定
  - _Requirements: 9.4, 9.5_

- [x] 5.3 マスク配列生成の実装
  - fill_mask!関数（各セルの内外判定、mask更新）
  - 複数物体への対応（任意数配置）
  - 物体内部速度の設定（デフォルト: [0,0,0]）
  - _Requirements: 9.1, 9.6, 9.7_

- [ ] 5.4 マスク生成の検証（マスク断面可視化）
  - テスト用Geometry JSONファイルの作成（box/cylinder/sphere各1個以上）
  - XY/XZ/YZ断面でのマスク値ヒートマップ描画
  - 物体形状が正しくマスク化されているか目視検証
  - 複数物体の重なり・隣接時の動作確認
  - _Requirements: 9.1, 9.4, 9.5, 9.6_

- [x] 6. 境界条件適用
- [x] 6.1 外部境界条件の実装
  - Dirichlet条件（境界値を直接代入）
  - Neumann条件（内部勾配をコピー）
  - 6面独立に速度境界条件を適用
  - 圧力のNeumann条件適用
  - _Requirements: 8.1, 8.2, 8.3_

- [x] 6.2 対流流出条件の実装
  - 輸送速度Ucの計算（境界面の平均速度）
  - 1次風上離散化（φ_new = φ - Uc·Δt·(φ - φ_neighbor)/Δx）
  - 各軸方向（x/y/z_min/max）への対応
  - apply_outflow!関数の実装
  - _Requirements: 8.6_

- [x] 6.3 吹出口・吸込口境界条件の実装
  - 矩形/円筒形状の吹出口・吸込口領域判定
  - 座標指定による速度設定
  - outflow/dirichlet条件の切り替え
  - _Requirements: 8.4_

- [x] 6.4 内部境界条件の実装
  - 部分境界（矩形/円筒領域）の領域判定
  - 法線方向に基づく境界面特定
  - マスク関数を用いた速度/圧力境界処理
  - _Requirements: 8.5_

- [x] 6.5 境界条件優先順位の実装
  - 優先順位ロジック（内部境界 > 物体 > 外部境界）
  - apply_boundary_conditions!関数の統合実装
  - _Requirements: 8.1, 8.2, 8.3, 8.4, 8.5_

- [ ] 6.6 境界条件適用のテスト
  - 単純格子での各境界条件タイプの動作検証
  - Dirichlet: 境界値が正しく設定されることを確認
  - Neumann: 勾配が正しく保存されることを確認
  - 優先順位: 重複領域で正しい条件が適用されることを確認
  - _Requirements: 8.1, 8.2, 8.3, 8.4, 8.5_

- [x] 7. 乱流モデル
- [x] 7.1 標準 Smagorinsky モデルの実装
  - 歪み速度テンソルの計算
  - 乱流粘性νt = (Cs·Δ)² · |S| の計算
  - Smagorinsky定数Csのパラメータ化（デフォルト0.2）
  - compute_turbulent_viscosity!関数の実装
  - _Requirements: 2.2, 2.3_

- [x] 7.2 有効粘性係数の計算
  - 分子粘性と乱流粘性の合算（nu_eff = nu + nu_t）
  - nu_eff配列の更新
  - _Requirements: 2.2_

- [x] 8. 物理モデル（対流・拡散）
- [x] 8.1 Lax-Friedrichsフラックス分割の実装
  - f±(u) = (1/2)(f(u) ± α·u) の計算
  - 各成分の最大速度α（max|u|, max|v|, max|w|）の計算
  - 各成分独立の再構成準備
  - _Requirements: 4.4_

- [x] 8.2 WENO3再構成（等間隔格子）の実装
  - Stencil 0, 1の候補多項式計算
  - 理想重み（d₀=1/3, d₁=2/3）の設定
  - 滑らかさ指標β_kの計算
  - 非線形重みω_kの計算（ε=10⁻⁶, p=2）
  - _Requirements: 4.2_

- [x] 8.3 WENO3再構成（非等間隔格子）の実装
  - 非等間隔格子再構成係数（a_k, b_k, c_k）の計算
  - 非等間隔理想重み（d₀, d₁）の計算
  - 代表セル幅Δx_refによる滑らかさ指標の無次元化
  - weno3_reconstruct_left関数の完成
  - _Requirements: 4.5_

- [x] 8.4 対流フラックス計算の実装
  - 3軸方向（x, y, z）の対流フラックス計算
  - f⁺（左から右）とf⁻（右から左）の再構成
  - 壁面境界でのマスク値による速度ゼロ処理
  - add_convection_flux!関数の実装
  - _Requirements: 4.1, 4.2_

- [x] 8.5 WENO3境界処理の実装
  - ゴーストセル設定（Dirichlet/Neumann/壁）
  - 境界付近での2次精度自動降格
  - _Requirements: 4.2_

- [ ] 8.6 対流項離散化のテスト
  - 等間隔格子でのWENO3再構成係数の標準値検証（d₀=1/3, d₁=2/3）
  - 滑らかな関数に対する再構成精度の検証
  - 不連続を含む関数に対する振動抑制の確認
  - 非等間隔格子での係数計算の正確性検証
  - _Requirements: 4.2, 4.4, 4.5_

- [x] 9. 拡散項離散化
- [x] 9.1 2次精度中心差分の実装
  - 非等間隔格子対応の中心差分離散化
  - 調和平均による界面粘性係数計算
  - _Requirements: 4.3_

- [x] 9.2 拡散フラックス計算の実装
  - 3軸方向の拡散フラックス計算
  - 壁面境界でのマスク関数による勾配ゼロ処理
  - add_diffusion_flux!関数の実装
  - _Requirements: 4.3_

- [ ] 9.3 拡散項離散化のテスト
  - 解析解のある熱伝導問題での精度検証
  - 壁面境界でのNeumann条件（勾配ゼロ）が正しく適用されることを確認
  - 非等間隔格子での拡散フラックス計算の正確性検証
  - _Requirements: 4.3_

- [x] 10. 圧力ポアソンソルバー
- [x] 10.1 Red-Black SOR法の実装
  - Red-Black分割によるSOR反復
  - 加速係数ωのパラメータ化
  - 相対残差ノルム ||b - Ax|| / ||b|| での収束判定
  - _Requirements: 7.1, 7.4, 7.5_

- [x] 10.2 (P) CG法の実装（オプション）
  - 共役勾配法の反復ループ
  - KrylovBuffers（r, p, q）の使用
  - 収束判定と反復回数管理
  - _Requirements: 7.2_

- [x] 10.3 (P) BiCGSTAB法の実装（オプション）
  - BiCGSTAB反復アルゴリズム
  - 追加作業配列（r0, s, t）の使用
  - 収束判定と反復回数管理
  - _Requirements: 7.3_

- [x] 10.4 圧力平均値引き戻しの実装
  - 流体セルの空間平均計算（物体セル除外）
  - 全流体セルから平均値を減算
  - solve_poisson!関数の完成
  - _Requirements: 7.6_

- [x] 10.5 収束失敗時の処理実装
  - DivergenceAction列挙型（WarnContinue, Abort）
  - PoissonConfig.on_divergenceに応じた動作分岐
  - 警告/エラーログ出力とチェックポイント保存
  - _Requirements: 7.4_

- [ ] 10.6 圧力ポアソンソルバーのテスト
  - 解析解のあるポアソン問題での収束性検証
  - 圧力平均値引き戻しの動作確認
  - 収束判定値の変更に対する反復回数の変化確認
  - 物体セル（mask=0）の除外が正しく機能することを確認
  - _Requirements: 7.1, 7.5, 7.6_

- [x] 11. Fractional Step法
- [x] 11.1 擬似速度計算の実装
  - 対流フラックス・拡散フラックスからの擬似速度計算
  - セルセンターでの(u*, v*, w*)計算
  - compute_pseudo_velocity関数の実装
  - _Requirements: 6.2_

- [x] 11.2 セルフェイス内挿の実装
  - 擬似速度のセルフェイスへの内挿
  - u*_{i+1/2} = (u*_i + u*_{i+1}) / 2
  - チェッカーボード不安定性の防止
  - interpolate_to_faces!関数の実装
  - _Requirements: 6.3_

- [x] 11.3 発散計算の実装
  - セルフェイス値を用いたセルの発散計算
  - div = (u*_{i+1/2} - u*_{i-1/2})/Δx + ...
  - ポアソン方程式のソース項生成
  - _Requirements: 6.4_

- [x] 11.4 速度補正の実装
  - 圧力勾配による速度補正
  - 連続の式を満たす速度場の生成
  - correct_velocity関数の実装
  - _Requirements: 6.5_

- [x] 11.5 Fractional Step統合の実装
  - fractional_step!関数の完成
  - 擬似速度→内挿→発散→ポアソン→補正のフロー統合
  - 境界条件適用タイミングの管理
  - _Requirements: 6.1_

- [ ] 11.6 Fractional Step法のテスト
  - 1タイムステップ実行後の速度発散が閾値以下であることを確認
  - セルフェイス内挿によるチェッカーボード抑制の確認
  - 速度補正後の連続の式の満足度検証
  - _Requirements: 6.1, 6.3, 6.4, 6.5_

- [x] 12. 時間積分スキーム
- [x] 12.1 Euler陽解法の実装
  - 1ステージ時間進行
  - euler_step!関数の実装
  - _Requirements: 5.1_

- [x] 12.2 (P) RK2法の実装
  - 2ステージRunge-Kutta法
  - RKBuffers（u_rk1, v_rk1, w_rk1）の使用
  - 中間値保存と最終更新
  - rk2_step!関数の実装
  - _Requirements: 5.2_

- [x] 12.3 (P) RK4法の実装
  - 4ステージ古典的Runge-Kutta法
  - RKBuffers全配列の使用（u_rk1〜u_rk4）
  - 勾配k1〜k4の計算と重み付け合成
  - rk4_step!関数の実装
  - _Requirements: 5.3_

- [x] 12.4 時間刻み計算の実装
  - CFL条件からの時間刻み自動計算
  - 拡散条件との整合性確認
  - compute_dt関数の実装
  - _Requirements: 11.9_

- [x] 12.5 時間積分統合の実装
  - advance!関数（スキーム選択に応じた分岐）
  - TimeConfig構造体による設定管理
  - _Requirements: 5.1, 5.2, 5.3_

- [x] 13. 出力機能
- [x] 13.1 SPHベクトル出力の実装
  - V-Isio SPH形式（svType=2, dType=1）のバイナリ出力
  - リトルエンディアン、単精度Float32
  - Fortran unformatted互換レコードマーカー
  - 有次元/無次元変換オプション
  - write_sph_vector関数の実装
  - _Requirements: 11.5_

- [x] 13.2 SPHスカラー出力の実装
  - 圧力SPH形式（svType=1, dType=1）のバイナリ出力
  - 瞬時値圧力ファイル出力
  - write_sph_scalar関数の実装
  - _Requirements: 11.5_

- [x] 13.3 ファイル命名規則の実装
  - vel_#######.sph, prs_#######.sph形式
  - vel_avg_#######.sph形式
  - 7桁ゼロ埋めステップ番号
  - generate_sph_filename関数の実装
  - _Requirements: 11.5, 11.13_

- [x] 13.4 チェックポイント出力の実装
  - 倍精度バイナリ、ヘッダ付き形式
  - ヘッダ（Nx,Ny,Nz, step, time, is_dimensional, L₀, U₀）
  - ゴーストセル込み配列データ
  - write_checkpoint関数の実装
  - _Requirements: 11.6_

- [x] 13.5 チェックポイント読込・リスタートの実装
  - チェックポイントファイル読込
  - 有次元データの無次元化処理
  - 代表量の一致検証（警告出力）
  - read_checkpoint関数の実装
  - _Requirements: 11.7, 11.8_

- [ ] 13.6 出力機能のテスト
  - SPHバイナリのエンディアン・レコードマーカー検証
  - 書き込んだSPHファイルの読み戻しによるデータ整合性確認
  - チェックポイント書込→読込の往復テスト
  - 有次元/無次元変換の正確性検証
  - _Requirements: 11.5, 11.6, 11.7, 11.8_

- [x] 14. モニタリング機能
- [x] 14.1 計算モニターデータ構造の実装
  - MonitorData構造体（step, time, dt, u_max, cfl, div_max, pressure_itr, pressure_residual）
  - MonitorConfig構造体（console_interval, history_interval, div_threshold）
  - _Requirements: 11.11_

- [x] 14.2 モニター値計算の実装
  - compute_u_max（速度最大値）
  - compute_cfl（実CFL数）
  - compute_divergence_max（最大発散）
  - _Requirements: 11.11_

- [x] 14.3 コンソール表示の実装
  - ステップ情報のフォーマット表示
  - 指定間隔での表示制御
  - _Requirements: 11.11_

- [x] 14.4 historyファイル出力の実装
  - テキスト形式（スペース区切り）
  - step time Umax divMax Pitr Presフォーマット
  - history.txtへの出力
  - _Requirements: 11.12_

- [x] 14.5 発散検出の実装
  - divMax閾値判定
  - NaN/Inf検出
  - エラー出力と計算停止処理
  - _Requirements: 11.16_

- [x] 15. 可視化機能
- [x] 15.1 可視化構成の実装
  - VizConfig構造体の定義
  - Package Extensions対応（CairoMakie weakdep）
  - 可視化無効時のstub関数
  - _Requirements: 12.1_

- [x] 15.2 断面抽出の実装
  - XY/XZ/YZ断面の抽出ロジック
  - plane_indexに基づく断面位置指定
  - 有次元量への変換
  - _Requirements: 12.3_

- [x] 15.3 コンター表示の実装
  - CairoMakieによるheatmap/contour描画
  - 速度場・圧力場のカラーマップ
  - カラーバー設定
  - _Requirements: 12.4_

- [ ] 15.4* (P) ベクトル表示の実装（オプション）
  - ベクトル矢印の重畳表示
  - 間引き数による密度調整
  - _Requirements: 12.5_

- [x] 15.5 画像出力の実装
  - PNG/SVG形式での画像保存
  - render_slice関数の完成
  - _Requirements: 12.6_

- [ ] 15.6* (P) テキスト出力の実装（オプション）
  - 断面データのテキストファイル出力
  - 座標+物理量のスペース区切り形式
  - write_slice_text関数の実装
  - _Requirements: 12.7_

- [x] 16. メインドライバ統合
- [x] 16.1 初期化フローの実装
  - パラメータ読込→格子生成→バッファ確保→マスク生成
  - 境界条件読込→初期条件設定
  - モニター初期化
  - _Requirements: 11.1, 11.2, 11.4_

- [x] 16.2 メインタイムループの実装
  - 時間積分ステップ実行（advance!）
  - 出力間隔判定と各種出力呼び出し
  - 時間平均更新（平均化開始時刻以降）
  - 発散チェックと異常終了処理
  - _Requirements: 5.1, 6.1, 11.11, 11.16_

- [x] 16.3 dry_runモードの実装
  - メモリ見積もり表示
  - 安定条件チェック（警告のみ）
  - 実計算をスキップして終了
  - _Requirements: 11.17, 11.18_

- [x] 16.4 リスタート機能の統合
  - start="restart"時のチェックポイント読込
  - ステップ数・時刻からの継続
  - _Requirements: 11.7, 11.8_

- [x] 17. 統合テスト
- [x] 17.1 単体テストの実装
  - WENO3再構成係数の等間隔格子標準値検証
  - 無次元化/有次元化の往復変換テスト
  - SPHバイナリのエンディアン・レコードマーカーテスト
  - 境界条件優先順位ロジックテスト
  - 圧力平均値引き戻しテスト
  - _Requirements: 4.2, 3.1, 11.5, 8.1, 7.6_

- [x] 17.2 統合テストの実装
  - JSON入力→格子生成→初期化フローテスト
  - Fractional Step 1タイムステップ実行テスト
  - チェックポイント書込→読込→リスタートテスト
  - 物体形状→マスク生成→境界条件適用テスト
  - _Requirements: 11.1, 6.1, 11.6, 11.7, 9.1_

- [ ] 17.3 E2Eテストの実装
  - 3Dキャビティフロー定常解検証
  - バックステップ流Re数依存性検証
  - _Requirements: 2.1_

- [ ] 17.4 性能テストの実装
  - 100³格子での1ステップ実行時間計測
  - メモリ使用量のスケーリング検証
  - _Requirements: NFR-4_

## Requirements Coverage

| 要件 | タスク |
|------|--------|
| 1.1-1.8 | 2.1, 2.2, 2.3, 2.4, 2.5 |
| 2.1-2.3 | 7.1, 7.2, 17.3 |
| 3.1-3.4 | 1.2, 2.4, 4.1, 17.1 |
| 4.1-4.5 | 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 9.1, 9.2, 9.3 |
| 5.1-5.3 | 12.1, 12.2, 12.3, 12.4, 12.5 |
| 6.1-6.5 | 11.1, 11.2, 11.3, 11.4, 11.5, 11.6 |
| 7.1-7.6 | 10.1, 10.2, 10.3, 10.4, 10.5, 10.6 |
| 8.1-8.6 | 6.1, 6.2, 6.3, 6.4, 6.5, 6.6 |
| 9.1-9.8 | 3.1, 5.1, 5.2, 5.3, 5.4 |
| 10.1 | 2.1 |
| 11.1-11.18 | 4.1, 4.2, 4.3, 4.4, 4.5, 13.1-13.6, 14.1-14.5, 3.5, 16.3 |
| 12.1-12.7 | 15.1, 15.2, 15.3, 15.4, 15.5, 15.6 |
| NFR-1, NFR-2, NFR-3, NFR-4, NFR-5 | 1.1, 2.1, 17.4 |
