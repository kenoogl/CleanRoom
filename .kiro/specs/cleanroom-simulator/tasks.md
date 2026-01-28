# 実装タスク

## 概要
クリーンルームシミュレータの実装タスク。既存実装と新仕様/設計の差分を解消し、全要件を満たす実装を完成させる。

---

## Task 1: 境界条件JSONスキーマの更新
**Requirements**: 8.1-8.8, 11.1, 11.2

InputReaderとBoundaryConditionsを新仕様スキーマに対応させる

- [x] 1.1: VelocityBCType列挙型を仕様準拠に更新（Wall/Symmetric/Periodic/Outflow/SlidingWall/Inflow/Opening）
- [x] 1.2: OpeningFlowType列挙型を追加（OpeningInlet/OpeningOutlet）
- [x] 1.3: Opening構造体を実装（name, type, boundary, position, size, center, radius, flow_type, velocity）
- [x] 1.4: InternalBoundary構造体を実装（物体扱い、速度は内部速度）
- [x] 1.5: BoundaryConditionSetを更新（openings配列でinlets/outletsを統合）
- [x] 1.6: InputReaderのJSON解析を新スキーマに対応（キーワード小文字正規化）
- [x] 1.7: 旧キーワード（dirichlet/neumann/inlet/outlet）を廃止、エラー出力
- [x] 1.8: Openingはopenings配列で扱い、ExternalBCでのOpening指定はエラーとする（優先順位: Opening > internal_boundaries > 外部境界）
- [x] 1.9: JSON全体のキーワード小文字正規化（start/dry_run/Time_Integration_Scheme/Poisson/on_divergence等、BC以外も対象）
- [x] 1.10: Inflow/Opening/SlidingWallの速度は一様固定（時間依存・分布プロファイル非対応）であることを明記・実装
- [x] 1.11: ExternalBCでOpeningを許可しない旨を設計書に明記（設計書反映）
- [x] 1.12: external_boundariesは`velocity.type/value`形式のみ許可し、文字列省略形はエラーとする

---

## Task 2: internal_boundariesの物体化
**Requirements**: 8.1-8.8, 9.1-9.8

internal_boundariesを物体（mask=0）として扱い、速度は物体内部速度として適用

- [x] 2.1: internal_boundaries読み込み時にGeometryと連携してmask=0を設定
- [x] 2.2: 物体内部速度としてvelocityを適用する処理を実装
- [x] 2.3: 境界条件優先順位の適用（Opening > internal_boundaries > 外部境界）

---

## Task 3: 圧力境界条件の仕様準拠
**Requirements**: 7.1-7.9, 8.2, 8.3

圧力境界条件をPeriodic以外は全てNeumann、Periodicは周期条件に統一

- [x] 3.1: apply_pressure_bcs!をNeumann/Periodic分離で再実装
- [x] 3.2: Outflow境界の圧力をDirichletからNeumannに変更
- [x] 3.3: 圧力平均引き戻しのコメント・処理をPeriodic例外対応に修正
- [x] 3.4: 圧力平均引き戻しは流体セル平均のみで行い、物体セルは平均値で埋める

---

## Task 4: 動的マスク処理の拡張
**Requirements**: 6.1-6.13, 8.1-8.8

圧力反復中にInflow/Outflow/Openingのマスクを一時的に0に変更

- [x] 4.1: マスク更新関数名を update_boundary_mask! に統一
- [x] 4.2: Inflow/Opening境界も動的マスク対象に追加
- [x] 4.3: FractionalStepで圧力計算前後にupdate_boundary_mask!を呼び出し
- [x] 4.4: 静的マスク方針を明記（Wall/Symmetric/SlidingWall=0、Periodic=1、Inflow/Outflow/Opening=1）

---

## Task 5: Outflow境界の取り扱い整理
**Requirements**: 4.1-4.9, 8.1-8.8

Outflow境界の粘性項は特別扱いせず、ゴーストセル値とマスクによって自然に評価する。対流流出の輸送速度は境界平均速度とし、擬似速度への適用時は **u^n** から算出する。

- [x] 5.1: Diffusionモジュールで境界条件情報を受け取る引数追加
- [x] 5.2: Outflow境界隣接セルの粘性項ゼロ化処理を削除（自然評価に統一）
- [x] 5.3: Outflow/Opening(outlet)の輸送速度Ucを **u^n** から算出して境界適用
- [x] 5.4: Outflow更新式の右辺は **u^n** を用いる（phi_refで評価）
- [x] 5.5: Outflow境界のセルフェイス擬似速度をセルセンター/セルフェイス値から構成（保存性重視）
- [x] 5.6: 圧力補正後にOutflow境界のセルフェイス速度を連続式から決定
- [x] 5.7: 初期条件でセルフェイス速度をセルセンター値から内挿して初期化

---

## Task 6: 弱圧縮性連続式の実装
**Requirements**: 2.1, 6.1-6.13

弱圧縮性を考慮した連続式（M² ∂p/∂t + ∂u_i/∂x_i = 0）を導入

- [x] 6.1: マッハ数M（= U_ref/c_sound）をパラメータとして追加
- [x] 6.2: ポアソン方程式のソース項にM²項を追加
- [x] 6.3: 圧力時間微分項の離散化を実装
- [x] 6.4: 弱圧縮性対応のためPoisson演算子をHelmholtz化（∇²p - αp）
- [x] 6.5: p^n保持用バッファを追加し、rhs構築で使用する
- [x] 6.6: Mach_number / Sound_Speed の入力仕様と検証を実装

---

## Task 7: ポアソンソルバーの仕様準拠
**Requirements**: 7.1-7.12

前処理sweep数とBiCGSTAB実装を仕様に合わせる

- [x] 7.1: PRECONDITIONER_SWEEPSを4に変更
- [x] 7.2: BiCGSTAB法を実装（前処理付き）
- [x] 7.3: on_divergence設定（WarnContinue/Abort）の動作を実装
- [x] 7.4: CG/BiCGSTABの前処理オプション（none/sor/rbsor/ssor）を追加し、noneは前処理なしとする
- [x] 7.5: CG/BiCGSTABをSPD形（A'=-A, b'=-b）で解く
- [x] 7.6: α=0の特異系でrhs平均値を除去・復元する
- [x] 7.7: 前処理後に周期境界を明示適用する

---

## Task 8: SPH/Checkpointフォーマットの仕様準拠
**Requirements**: 11.1-11.18

SPH出力とチェックポイントのフォーマットを仕様に合わせる

- [x] 8.1: SPHのDZを非等間隔時Lz/Nzに変更
- [x] 8.2: SPHベクトル出力をu,v,wブロック順に変更（インターリーブ廃止）
- [x] 8.3: Checkpointのエンディアンをリトルエンディアン固定に変更
- [x] 8.4: is_dimensional検証（0以外はエラー）を追加
- [x] 8.5: Reference_Length/Velocityのリスタート時一致確認を実装（不一致時は警告）
- [x] 8.6: SPHは常に有次元量で出力（velocity/pressure/vel_avg）
- [x] 8.7: SPHレコード構成（Time/Step/Origin/Pitch）を仕様通りに固定
- [x] 8.8: SPH出力を単精度（Float32）・Fortran unformatted互換で確認

---

## Task 9: 時間平均機能の実装
**Requirements**: 11.13, 11.14

Start_time_for_averaging以降のWelford法による時間平均と出力

- [x] 9.1: update_time_average!をメインループから呼び出し
- [x] 9.2: 平均化開始時刻の判定処理を追加
- [x] 9.3: vel_avg SPHファイル出力を実装

---

## Task 10: 可視化の外部化
**Requirements**: 12.1-12.8

内蔵可視化を廃止し、外部ツールに移行する

- [x] 10.1: 内蔵可視化コード（Visualization.jl / Extension / 呼び出し）を削除
- [x] 10.2: ソルバーJSONからVisualizationブロックを削除し、存在時は警告する
- [x] 10.3: 可視化設定JSONを用いた実行経路（dispatcher / 各ツール）を追加
- [x] 10.4: 旧内蔵断面可視化を `visualize_cavity.jl` の slice モードへ移植
- [x] 10.5: ツールREADMEを更新し、JSON設定例を追加
- [x] 10.6: 可視化出力先をJSONのoutput_dirで制御（options.output_dir優先）
- [x] 10.7: cavityのslice実行時に中心線プロファイルをoutput_dirへ出力

---

## Task 11: 毎ステップ安定条件監視
**Requirements**: 5.1-5.3, 11.9, 11.10

計算中も各ステップでCFL/拡散数を監視し、条件違反時は停止

- [x] 11.1: compute_cfl_max/compute_diffusion_maxをメインループで呼び出し
- [x] 11.2: 条件違反時の警告・停止処理を実装
- [x] 11.3: Monitorにcfl_max出力を追加
- [x] 11.4: dry_run時は警告のみで停止しない挙動を明記・実装
- [x] 11.5: debug有効時にInflow/Outflow境界面の流量モニタ（正味流量、out/in分離、un最小/最大）を出力

---

## Task 12: WENO3右側再構成の実装
**Requirements**: 4.1-4.9

f⁻用の右側再構成関数を追加

- [x] 12.1: weno3_reconstruct_right関数を実装（左右反転ステンシル）
- [x] 12.2: add_convection_flux!で右側再構成を使用
- [x] 12.3: WENO3の入力がセル平均であることを明記・確認
- [x] 12.4: 非等間隔格子対応のWENO係数を実装・検証

---

## Task 13: condition.txt出力の実装
**Requirements**: 11.1-11.18

シミュレーション開始前に計算条件を出力

- [x] 13.1: condition.txt出力関数を実装
- [x] 13.2: 物理パラメータ、格子、時間積分、境界条件、出力間隔を記録
- [x] 13.3: condition.txtにOpening/内部境界/境界優先順位/動的マスク設定を含める
- [x] 13.4: condition.txtに前処理オプション（none/sor/rbsor/ssor）を出力する

---

## Task 14: 統合テストとバリデーション
**Requirements**: 全要件

全機能の動作確認とベンチマーク検証

- [x] 14.1: 新スキーマJSONでの読み込みテスト
- [x] 14.2: 3Dキャビティフロー検証ケースの実行
- [x] 14.3: SPH/Checkpoint出力の検証
- [x] 14.4: リスタート機能の検証
- [ ] 14.5: ポアソン係数行列の対称性を確認（数値検証またはコードレビュー）

---

## Task 15: Z方向非等間隔格子対応
**Requirements**: 1.3, 1.4, 1.5

Z_grid.type=non-uniform時のZ座標ファイル読み込みと格子生成

- [x] 15.1: read_z_grid_file関数の実装（格子点数=Nz+5の検証含む）
- [x] 15.2: Z_grid.type切り替えロジックの実装（uniform/non-uniform）
- [x] 15.3: 非等間隔格子でのWENO3係数計算の検証テスト
- [x] 15.4: SmagorinskyのΔを局所セル幅で評価する実装と検証
- [x] 15.5: non-uniform時のLz算出（z_{Nz+3}-z_{3}）とDZ=Lz/Nzの保持

---

## Task 16: SlidingWall境界条件の実装
**Requirements**: 8.1

壁面スライド駆動（接線速度指定、法線速度ゼロ）の実装

- [x] 16.1: 境界面ごとの速度成分解釈を実装（法線成分無視、接線成分適用）
- [x] 16.2: ゴーストセルマスク=0の設定
- [x] 16.3: 壁面せん断流束の補正をDiffusionモジュールに追加

---

## Task 17: RK2/RK4時間積分スキーム
**Requirements**: 5.2, 5.3

Euler以外の時間積分スキームの実装

- [x] 17.1: RKBuffers構造体の条件付き確保ロジック
- [x] 17.2: rk2_step!（Heun法）の実装（Butcher tableau準拠）
- [x] 17.3: rk4_step!（古典的4段RK）の実装
- [x] 17.4: Time_Integration_Scheme切り替えロジックの実装

---

## Task 18: dry_runモード
**Requirements**: 11.17, 11.18

ドライラン時のメモリチェックと実行可否判定

- [x] 18.1: estimate_memory_size関数の実装（CFDBuffers + Krylov + RK配列）
- [x] 18.2: Sys.free_memory()/Sys.total_memory()との比較
- [x] 18.3: 警告出力と実行可否判定の表示

---

## Task 19: history.txt出力拡張
**Requirements**: 11.11, 11.12

divMax位置(i,j,k)の出力とフォーマット調整

- [x] 19.1: compute_divergence_max関数に位置情報の返却を追加
- [x] 19.2: MonitorDataにdiv_max_pos::NTuple{3,Int}を追加
- [x] 19.3: history.txt/コンソール出力フォーマットの更新

---

## Task 20: 逆流安定化オプション
**Requirements**: 8.1

Outflow境界での逆流時にフラックスをゼロクリップ

- [x] 20.1: 逆流検出ロジックの実装（法線方向速度の符号判定）
- [x] 20.2: 対流フラックスのゼロクリップ処理
- [x] 20.3: オプションパラメータ（reverse_flow_stabilization）の追加

---

## Task 21: 可視化ツール整理（問題クラス別パターン）
**Requirements**: 12.1-12.8

問題クラスごとの可視化ツールを整理し、適用パターンを明確化

- [x] 21.1: cavity/backward_step のツールをパターンとして整理
- [x] 21.2: パターン適用方法（CLI/引数）をドキュメント化
- [x] 21.3: 統一入口（dispatcher）を追加し、パターン追加方針を明記

---

## Task 22: RBSSOR（rbsorによる4スイープ対称構造）の実装
**Requirements**: 7.1

RBSSORソルバーとPrecondRBSSOR前処理の実装（4スイープ対称構造）

- [x] 22.1: SolverType/PreconditionerTypeにRBSSOR/PrecondRBSSORを追加
- [x] 22.2: sor_sweep_rbssor!関数を実装（前進R→B, 後退B→R, 前進B→R, 後退R→B）
- [x] 22.3: solve_poisson!/apply_preconditioner!でRBSSORをディスパッチ
- [x] 22.4: InputReader/MonitorにRBSSOR対応を追加
- [x] 22.5: 設計書・要件書にRBSSORを追加
- [x] 22.6: RBSSOR収束性検証（RBSOR: 40反復, SSOR: 35反復, RBSSOR: 53反復で収束確認）

---

## Task 23: 前処理の符号整合性（negate_rhsフラグ）

**Status**: COMPLETED

**Requirements**: 7.5

CG/BiCGSTAB の SPD 形式と SOR 前処理の符号整合性を確保する

- [x] 23.1: sor_sweep_rbsor! に negate_rhs パラメータを追加し、b_val = rhs_sign * rhs * vol で計算
- [x] 23.2: sor_sweep_lex! に negate_rhs パラメータを追加
- [x] 23.3: sor_sweep_rbssor! に negate_rhs パラメータを追加
- [x] 23.4: apply_preconditioner! で negate_rhs=true を指定して呼び出す
- [x] 23.5: CG + SSOR 収束検証（以前発散→27反復で収束）
- [x] 23.6: BiCGSTAB + SSOR 収束検証（17反復で収束、継続動作確認）
- [x] 23.7: 設計書に符号整合性の説明を追加

---

## 優先順位（推奨実行順序）
1. Task 1: 境界条件JSONスキーマ（最重要、他タスクの基盤）
2. Task 2: internal_boundariesの物体化
3. Task 3: 圧力境界条件
4. Task 4: 動的マスク処理
5. Task 7: ポアソンソルバー
6. Task 5: Outflow境界の取り扱い整理
7. Task 6: 弱圧縮性連続式
8. Task 8: SPH/Checkpointフォーマット
9. Task 12: WENO3右側再構成
10. Task 9: 時間平均
11. Task 10: 可視化（外部化）
12. Task 21: 可視化ツール整理
13. Task 11: 安定条件監視
14. Task 13: condition.txt
15. Task 14: 統合テスト
16. Task 15: Z方向非等間隔格子
17. Task 16: SlidingWall境界条件
18. Task 17: RK2/RK4時間積分
19. Task 18: dry_runモード
20. Task 19: history.txt出力拡張
21. Task 20: 逆流安定化オプション
22. Task 22: ssor前処理のrbsor実装
23. Task 23: 前処理の符号整合性
