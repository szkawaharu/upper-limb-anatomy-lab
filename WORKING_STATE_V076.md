# v0.76 作業状態（2026-11-27、継続中）

公開mainはc81256f / v0.75.0。作業ブランチredesign/v076-shoulder。まだ公開しない。ユーザーはv0.76を完成し実装することを依頼している。骨格修正は検証が進んだが筋の品質が未達。mainへ無理に出さないという元の条件を維持する。

## 確定した修正

- 実肋骨を包む162支持面、source肩甲骨987点、文献priorを基準にした接平面の連成拘束。反復中に包絡面を切り替えない。48反復＋バックトラック。
- 不均一C1関節補間1,408節点、11,907検証点、最大剛体点誤差上限0.967092mm。誤差1mmを変更していない。参照依存SHA256も記録。
- 543姿勢（0/40/90°面、挙上0〜180°各1°）で56,630元肩甲骨頂点の肋骨内部0。代表21姿勢の肩甲骨辺→肋骨面／逆方向の辺→面の交差0。共面重なりと他関節は未評価。
- 旧13姿勢の骨格依存ウェイト選択から共通基準フィールドを分離。注意：v61SelectWeightsは4引数のshoulderPoseFor(..., true)を呼ぶため、sync=trueでありskipGuide=trueではない。旧記録の骨格と基準データの結合の指摘は有効。
- 最大辺ペア・最大面積の三角形IDと元の物理寸法を監査へ追加。診断APIは操作キューをflushし最終骨姿勢を読む。
- 32試験通過後、監査単位の追加試験と診断API修正を入れ、33試験を再実行して全通過。48対96反復の184姿勢比較も最大0.0000142mmで通過。ログtest-results/v076-source-bone-offline-33.log。

## 筋の根本原因・データ

- 元STLの閉鎖性はクラスタリングで壊れ、重複逆向き面の片方だけを捨てる処理も穴を作っていた。exact weld、逆向き対を両方相殺、未参照頂点の除去で外形を変更せず修復できる。sourceVertexIdsを保持。
- 菱形筋の0/1固定領域が両方欠如、多数区分で停止側が欠如。旧ランドマークとsource筋／骨表面の位置にずれがある。docs/validation/v076/source-bones/legacy-attachments.json。
- scripts/calibration/infer-terminal-regions.pyは実source骨表面と終端候補を照合。骨の実測エッジ間隔を距離許容に使う。肩甲骨の内側縁ROIと後面TS→AA経路を利用。遠い点を選んだ局所ビン案は取り消し、現行は近い候補だけ。まだ臨床・手動解剖照合に合格したデータではない。
- 原型source全点のデータは以下（大きいためtest-results内、Git追跡外）。いずれも同じworld座標、scale=.00555/mm。
  - v076-material-baseline/bindings.json：57e時点の旧クラスタ基準、約70MB。
  - v076-source-bindings/bindings.json：exact source weld、ただし旧の逆向き面片方除去の欠陥あり。
  - v076-source-clean/bindings.json：面修復後、元頂点対応付き、約65MB。
  - v076-source-terminals/bindings.json：欠けた付着候補を追加したもの。PecAbdの起始は腱膜sourceが未照合で欠如。既存領域保持。
  - v076-source-attached/bindings.json：上記のDirichlet調和ウェイト。両端のある成分のみCGを解き、他成分は旧値保持。source固定マスク一致。
  - v076-source-refit/bindings.json：新調和場に整合する正の重みのcentroid Bézier control polygonsへrefit。material=tは調和場を正規化。別投影後平滑化を廃止する実験。
  - v076-source-refit/material/：静的中心、source法線レイ半弦ベクトルの80回平滑場、体積積分の層別求積。
- 厚みレイは三角形端面で長い値となる。実測筋厚さと呼ばない。新規形状を捏造しない。

## 不採用・未統合の変形実験

scripts/experimental以下は製品既定に統合していない。
- 単独closest material projection：全体の辺・面積外れ値悪化（前回履歴）。
- 旧ガイドの直線化、中心回転、DQSだけ：外れ値・体積の問題を解消せず。
- rest-surface.js／rest-surface-v1-rejected.js：局所微分面積＋半弦厚み補正。全体悪化。不採用。
- rest-transport.js／surface-transport-volume-v1.js：中心回転またはchord transport＋閉鎖source体積のスカラー投影。半弦ベクトルはsource法線に沿う計測を平滑化し、手動筋別係数はない。大きい体積喪失・負の符号付き体積が残り、投影の根が求まらない姿勢がある。制限を緩めて通していない。
- audit-refit-engine.mjs：新付着／調和場／正のガイドに旧共通エンジンを適用。中立・90°は改善するが高挙上で体積0.24〜0.55などになり不合格。source-refit/common-engine-audit.json。
- audit-volume-transport.mjs：source-refit/material/transport-audit-{chord,surface}.json。単純な回転補間と全体体積だけでは局所折返しを直せない。
- tests/helpers/source-preview.mjsはテスト専用CPU描画のhydration。QA_PREVIEW_SOURCE=1でのみ使用。これを製品の合格済みエンジンと誤認しない。
- 比較画像test-results/v076-source-center-preview/：57eとsource-centerの同条件18画像。見た目が滑らかでも体積不合格なので不採用。

## 次の実装判断

正しい固定付着を境界条件とし、局所フレームの特異点を持たない共通の変形を検討・実装する。単に筋厚みや許容値を変えない。変形グラフ／ARAP等の幾何学的な共通拘束を調べるところ。まだグラフエンジンは実装していない。実装時はsource全解剖メッシュを残し、coarse graphが形状を置換しないようにする。骨最終姿勢・固定領域・左右／首の独立性をテスト。体積・接触・局所破綻を全体で測定し、最後に同条件の画像比較・main／Pages更新を行う。改善が未達なら完成・公開済みと偽らない。

## 実行環境・再現

- Node/Playwright、Chrome /Applications/Google Chrome.app/Contents/MacOS/Google Chrome。
- source STLキャッシュ /tmp/upper-limb-atlas。Python /tmp/upper-limb-research-venv/bin/python（numpy/scipy/matplotlib/trimesh/rtree/embreex）。システムpython3にはnumpyがない。
- 全ブラウザQAは専用一時HTTPサーバー。実機スマホ未検証。
- npm test / node scripts/calibration/check-shoulder-reference.mjs / npm run build:pages。
- node scripts/calibration/build-shoulder-field.mjs tests/fixtures/shoulder-v076-field-calibration.json：test-results/v076-source-envelope-reference-cache.jsonを依存SHA付きで再利用。初回は時間がかかる。
- QA_NUMERIC_ONLY=1 QA_EXPORT_BINDINGS=1 QA_SOURCE_MESH=1 QA_OUTPUT=... node tests/browser-v073.mjs：元source exact weldの診断export。
- scripts/calibration/restore-source-topology.py → infer-terminal-regions.py → build-harmonic-attachments.py → refit-source-guides.mjs の順。
- `QA_BONE_FIELD=1 QA_BONE_SWEEP=1 node scripts/calibration/audit-bone-surfaces.mjs outdir`、Python audit-bone-surfaces.pyで実骨内外。QA_BONE_EDGES=1で双方の辺→面も検査。

GitHubは利用可能。mainへforce reset/pushしない。公開build-info SHAは最後にfetchして確認。ユーザーから実装・GitHub・公開は許可済みで、通常判断の確認は不要。サブエージェントは開発者指示により禁止。
