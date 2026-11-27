# v0.76 作業状態（2026-11-27、継続中）

公開mainはc81256f / v0.75.0。作業ブランチredesign/v076-shoulder。まだ公開しない。ユーザーはv0.76を完成し実装することを依頼している。骨格修正は検証が進んだが筋の品質が未達。mainへ無理に出さないという元の条件を維持する。

## さらに新しい状態：実距離・実関節列・材料積分

最新版は[HANDOFF](HANDOFF_V076.md)、[JOINT_CONTACT_V076](docs/JOINT_CONTACT_V076.md)。旧求積の実関節列は90→147.375°、55収束姿勢、固定・非正J・面向き警告0だが表面隣接J=0.014495〜10.489311で停止。180°と下降は未到達。全体積比約1だけで採用しない。

新しい凸交差積分は元51,352四面体と18,124制御セルを照合。被覆相対差1.782e−12、最大元要素相対差0.01361%、交差量0の元要素0。背景の旧正則化密度・元表面・固定面・骨は維持。90°で同じ実距離目的関数の対照と再収束を比較し、表面隣接Jは0.107479〜7.572327→0.252048〜3.273161。改善しても局所体積は未合格。新候補はtest-results/v076-mid-intersection-integrated-v3、対照はv076-mid-original-integration-control。

120°で12同条件画像を取得し両パネルを目視、例外／GLエラー0。他の筋は旧方式で不良が残る。数値・画像はdocs/validation/v076/joint-contact。Pythonは追加3件を含む23試験。新積分の高挙上・順逆・全筋・本番CPU/GPU・実機は未実施。下記は前段階の履歴。

## 最新状態：共通接触・原型の連続性参照

詳細と再開手順は[HANDOFF_V076.md](HANDOFF_V076.md)の先頭、[CONTACT_REFERENCE_V076.md](docs/CONTACT_REFERENCE_V076.md)。この文書の下方は時系列履歴で、古い「次」は実施済みの場合がある。

実骨接触を共有材料ソルバへ追加し、三角筋中部90°で新規深さ増加を2.281→0.043mm、全232,262表面サンプルでは最大0.057mmへ減らした。全体積0.999107、固定0。ただし元からの約4.075mmの骨内重なりは残り、全可動域と自己接触は未合格。

自己交差のある僧帽筋上部の別修復参照を作成。実際の3°/6°骨姿勢で旧点168と同じ空間位置の移動4.977mm、0.375°では0.704mm。ただし追加停止候補・修復原型を含み、製品の旧曲線19.1166mmが直ったとはしない。bulk3000でも6°のJは0.870〜1.162で局所体積は未達。

ローカル36 Node＋20 Python代数＋実原型回帰を通過。後者は正J・指定骨の固定・全体積・反復中立復帰（ドリフト1.524e−11mm）を確認。追加CIの実行結果はコミット後に照合。画像24枚とパネル、全数値はdocs/validation/v076/contact-reference。新しいWeb変形エンジンの統合、全筋・順逆・高挙上・左右・頸部・GPU・スマホ30fpsは未完了。

次は実際の関節角列で接触と正Jを同時に維持する継続。接触なしの180°途中解progress0.703733を完成姿勢と誤用しない。ユーザーの継続指示は有効で、品質未達を許可待ちや環境不足と説明しない。

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

正しい固定付着を境界条件とし、局所フレームの特異点を持たない共通変形を実装する。単に筋厚みや許容値を変えない。現在はARAP変形グラフの参照実装まで作成したが、以下の理由で不採用。source全解剖メッシュをグラフへ置換していない。全体体積を保つだけでは局所折返しを直せない。固定面全体との連続な接続、局所の非反転、骨面の滑走を同時に扱う必要がある。

## 最新の追加作業（2026-11-27、f71a78f以降）

- 共通ARAPはscripts/experimental/deformation-graph.js。局所回転の対称固有値解、固定境界、正の辺エネルギー、事前Cholesky、4節点の埋込。約200〜344節点。全source積分による全体体積の拘束も実装したが、面の向きの警告が残るため製品へロードしない。
- 三角筋中部の代表6姿勢では全体体積比はほぼ1でも、180°外転で274面、180°前方で99面に移送法線に対する面方向警告が残る。体積ヤコビアンや厳密自己交差は未評価。反復内の全source積分は秒単位になる場合があり、そのままではリアルタイム用に採用できない。
- 固定領域の隣の変形を四次基底で接続する案、DQS初期値、回転を正規化した埋込も全体を合格にしなかった。DQS単独＋体積投影も負の符号付き体積と根がない姿勢が残る。境界係数・厚み・許容値を緩めて通していない。
- グラフ生成の大胸筋腹部には低信頼の下端cap起始仮説がある。腱膜のsource照合を終えたデータではなく、製品へコピーしない。各筋の特性は定義データへ分離してから採用判断する。
- 合成テストで既知の剛体回転・ほぼ平面の回転・固定領域・全体積の微分・中立復帰・保存済み左右姿勢の非干渉を確認。参照状態が次の骨計算に書き換わる共有バッファ問題も防止した。医学的合格を意味しない。
- 全34オフライン試験を通過。f71a78fのGitHub CI（45850937304）も成功。作業ブランチもCI対象に追加した。Pagesジョブはmain以外ではskipする仕様。
- 実アプリの最新監査はdocs/validation/v076/source-bones/browser-numeric/。3°順逆掃引の最大サンプル移動19.116553mmは僧帽筋上部168、外転3→6°。0.375°へ細分すると最大3.018282mm（比0.1579）。大きい連続変形が残る。優先6区分の7姿勢では最大辺長比18.19877、最大面積比1335.02460（三角筋後部、前方90°）。旧21.165mmを最新最大と呼ばない。
- 同ブラウザ監査の復帰誤差0、検査対象NaN0、ブラウザエラー0、GPU18,413点のCPU差最大0.00277081mm。頭部ドラッグ・左右独立性のassertも通過。形状・体積の合格とは区別する。
- 公開v0.75 c81256fと作業版の同カメラ／同UI角で36画像を撮影し目視確認。120°外転、120°前方、両腕180°＋首30°、骨／筋、正面／側面／背面。docs/validation/v076/source-bones/comparison/に原画と6パネル。筋の肩周囲の形状不良が残る。
- 最終SC/AC/ST角とprior差をsource-bones/current-motion.jsonへ保存。180°のST上方回旋は54〜57°だが、priorから内旋約14〜18°増、後傾約12〜14°減の差もある。元loft項だけの比較も保存。実骨の接触0だけで解剖学的合格としない。
- 保存済み基準データの個別診断では肩甲挙筋78の前方90→93°は2.13463mm、僧帽筋下部6696の外転120→123°は1.31049mm。現在の全筋最大値とは別の固定入力比較。frame-traces/にrest／posedフレームと骨姿勢を保存。現在の僧帽筋上部168は、基準の横断幅が1mm未満の悪条件チャートから大きく増幅される。alignの反対向き切替には該当しない（dot最小0.9968）。

新しい採否の説明はdocs/MUSCLE_DEFORMATION_V076.md、再生成コマンドはscripts/experimental/README.md。大きいbindings・グラフ・f32はtest-results内でGit追跡外。採否に必要な集計・反例・入力SHAはdocs/validation/v076/muscle-experimentsへ保存した。次の通常チャットは、未採用の試作を既定へ昇格させたり、過去の数値を現在の最大として引用したりしない。

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

## 次の継続で追加した試作（2026-11-27朝）

- チャートの全2方向微分とsource辺間隔によるridge逆写像を試した（differential-chart.js）。僧帽筋上部168の保存例は19.1166→7.4351mmへ減るが、元表面の大胸筋・三角筋全体は悪化し、体積比0.135や辺長比92超が残るため不採用。原型・旧クラスタの両方を比較し、単一点の改善を採用理由にしていない。
- volumetric_reference.pyは曲線を使わない四面体の静的参照。三角筋中部の全原型頂点・面を保持する51,352四面体を生成し、source境界完全一致・正体積・全体積一致を確認。合成四面体の独立数値微分・剛体不変性・反転拒否テスト通過。
- 初回基本肢位：体積比1.000450、非正要素0、付着誤差0、辺長比0.767〜1.264。ただし局所J=0.377〜2.184が残り、停止残差もbulkでスケールしていたので強化が必要。90°試行はその見直しで途中停止。未完了を合格としない。docs/VOLUMETRIC_REFERENCE_V076.md。
- 大きい参照データはtest-results/v076-volumetric。Python venvへ研究ツールtetgen0.8.4とPDF閲覧用pypdfium2を追加したがWebランタイム依存にはしていない。TetGen wrapper MIT／core AGPLv3を記録。生成物の元アトラスの出典を維持し、バイナリを公開物へ入れない。
- Liu2017の著者PDFは/tmp/liu-2017-quasi-newton.pdf、本文は同.txt。p5〜7のアルゴリズムとp10/12を確認。動的シミュレーションや実測材料定数を再現したとは主張しない。

## 体積参照の継続（2026-11-27、cbd03a0以降）

- 停止判定をbulk／L-BFGSスケールに依存しない `L_FEM^-1 g < 0.001mm` へ変更。近似SPDの行列不要準Newton法と、独立した微分・密行列・回転共変性テストを追加した。
- 元51,352要素の基本肢位をbulk=3000で再収束：全体積比1.000005197、局所J=0.796296〜1.666430、非正要素0、固定端変化0、残差0.00080796mm。bulk=30でも収束したがJ=0.3764〜2.1851。収束後も局所非圧縮性の外れ値が残る。
- 元の全境界頂点・面を保持した別の内部分割（16,006点、56,083要素）を生成。境界完全一致・体積一致。条件数中央値21.01→12.02、最大値は14819.61のまま。旧解を内点へ移しただけでは6要素が非正になったため転送を拒否した。新分割は固有の境界継続から解き直した。
- 新分割の基本肢位bulk=3000は30外側反復でJ=0.7550〜1.3191、残差0.002199mmで未収束。小さい体積外れ値と反復不足を混同しない。スカラー前処理で追加してもCG上限250へ到達する。source-rest連成前処理も試作したが、最初の一反復で上限へ達している。角度や閾値は変更していない。
- 24制御点の調和アフィンパッチ96スカラーモードは、進捗0.54744でJ最小0.009215、全FEM残差23.47mm。停止した不採用の途中結果。表面スカラーウェイトだけのliftingは中立のcold startで80負要素を作った。4チャンネルの体積FEM調和liftingへ変更しても、少ないモードの変形能力は不足する。
- sourceの厳密な線形弾性接線応答から12ベクトルモードを作成。独立微分、剛体の無ひずみ、固定端0、密行列、完全空間への射影テストは通過。しかし実メッシュの中立へ向かう進捗0.749707でJ最小0.02676、全FEM残差47.04mmとなり停止。不採用。小さな射影勾配や全体体積約1だけを採用理由にしない。
- 新しいPythonの5試験をCIへ追加（numpy/scipyのみ、TetGenはCI・Web依存に含めない）。全35 Node試験と区別する。CI結果は次のコミットのActionsを確認すること。
- 大きいデータは `test-results/v076-volumetric` と `test-results/v076-volumetric-quality`。研究参照コードはscripts/experimental。Webはこの試作をロードしていない。筋・接触・性能の合格は未達でありmainへ出さない。

## 収束比較・新しい引継ぎ（b08e40f以降）

- b08e40f57aa7b490366d51233840b98d0f0e1af8をpush。Actions 45907316836で35 Node＋5 Python、Pages組立検査が成功。作業ブランチの公開ジョブはskip。mainと公開はc81256fのまま。
- source-rest連成前処理、節点polar回転で運ぶ前処理はいずれも実データでCG250上限へ達し、途中停止した。現在の変形勾配の連成SPD行列を組んで因子化する `frozen` 前処理へ変更。独立した行列不要演算と一致。これは線形解法の変更だけで、エネルギー・固定端・停止条件を変えない。
- 新分割基本肢位を厳密残差0.000599213mmへ再収束。J=0.7550179〜1.3191062、体積比1.0000087638、固定端変化0、非正要素0、Float32後も非正0。辺長比0.322990〜1.693932、面積比0.264614〜1.984910。局所形状・体積はまだ未達。
- 元分割と新分割の全source表面差は最大0.659523mm／RMS0.354994mm。固定領域差は丸め誤差のみ。最大Jの要素の条件数は約81.58で、最高条件数sliverだけに起因する問題ではない。分割2例によるメッシュ収束性の証明としない。
- 現在の収束解：`test-results/v076-volumetric-quality/neutral-refined-gn-k3000-frozen-v2.npz/json`。旧分割は`test-results/v076-volumetric/neutral-refined-gn-k3000.npz/json`。compact結果と独立監査はdocs/validation/v076/muscle-experiments/continuum-v2。
- 同じsourceと骨姿勢の静止比較をtests/browser-v076-solid.mjsへ追加。最初の旧分割previewは重複した体積集計のTypedArray.map誤りで追加欄がnull（共通監査は有効）になったため採用しない。修正済みコードと新分割で再撮影する。`test-results/v076-solid-neutral`を正式な新比較と誤用しない。
- HANDOFF_V076.mdを現状中心に書き換えた。古い21.165mmを最新最大とする依頼文や、実骨未監査とする再開手順は廃止。継続時はこの指示書と既知の問題を読む。

- 修正後の新分割比較は完了：6組12画像、最終骨a/b・カメラ一致、例外0／GLエラー0。全筋／診断用の三角筋中部、正面／側面／背面。両パネルを目視確認。輪郭の差は小さく、挙上の改善とは扱わない。曲線方式に比べ体積比と最大辺／面積は下がるが最小辺／面積は悪化する。docs/validation/v076/solid-neutralに入力SHAと原画像を保存。描画runtimeがb08e40fから不変なことも照合した。
- 状態依存連成行列の追加後、5 Python試験を再実行し全通過。独立の4要素共有節点組立、現在変形の行列不要演算との一致も追加確認。35 Nodeはこの継続で再実行し全通過。以降の変更は研究コード・ブラウザ比較試験・文書だけであり、製品runtimeは変更していない。

## 2026-11-27：今回の続行（fb31d24以降、まだ作業中）

### 続行の最新（da4182f以降、詳細はdocs/ATTACHMENT_BOUNDARIES_V076.md）

- da4182fe7ac395efd2223178b8b7c544c8155677を作業ブランチへpush、CI45939904659で35 Node＋10 Pythonが成功。main／公開はその後もc81256f。ユーザーの許可不足ではなく品質未達で未公開。
- 実骨の肩甲骨→上腕骨内判定は90°外転・前方とも全56,630点で0。逆向き・全面交差の合格とは区別。筋の旧固定点が骨内に入る矛盾は残る。
- 外部モデルは骨ローカルの相似登録を試し、TS/AA/AI残差7.8〜12.3mm、参照AAの骨面距離15.38mmを確認。実骨ランドマークと同一とみなして座標移植しない。source model geometryのローカルcacheは/tmp/v076-reference-model。登録データはtest-results/v076-reference-registration。
- 起始候補を肩峰外側壁の安静時ROIと実メッシュ解像度から497点へ再分類（旧392点と共通57点）。実測enthesisではない。骨面距離0.00049〜1.54715mm、帯幅は実辺長から1.55108mm。全元表面と停止推定496点を保持。候補＋粗18mm格子の90°は体積0.997517、面向き警告0、肩甲骨への新規交差0だが上腕骨の新規深さ4.00902mm。
- **共通エンジンの新原因**：粗格子の支持点固定が未宣言の組織4,313点まで固定する。最悪新点28359は重み0.0975でも支持3点が全て起始固定（1点の自由支持は重み0）。基準点12473も99.6%固定。単に付着マスクを改善しても、補間固定範囲が膨張する。
- 通常のDelaunayに元固定993頂点を入れるだけでも、付着面内が自由制御点で補間される反例がある。実データ検査で拒否。point_cage.pyは未採用として保存、合成の点拘束試験だけを合格理由にしない。
- PLCで元993固定点・1,702固定面を回復。TetGen AGPLはオフライン生成のみ。未宣言剛体化は57新表示点に減る。90°は全段階収束413.6秒、体積0.999074、面向き0、固定誤差0、肩甲骨への新規0、上腕骨の新規最大2.27577mm。接触なしでは不合格。
- attachment_refinement.pyは安静時の宣言外固定の制御辺／面を細分。源の条件から17辺、80面、84セルに計181の自由制御点を追加。source-boundary-directed refinement後3,044点／18,124要素、38,712表示点の宣言外剛体化0。元表面・付着面・親セル積分密度を保持。全域細分は11,332点／54,825要素になるので標準にはしない。選択は接触・姿勢・見た目に依存せず、材料や解剖の係数を足していない。
- 最新候補test-results/v076-plc-surface-refined-mid-acromial-18の90°解き直しを実行中。ソルバは厳密global Newtonを下降方向・正J・Armijo付きで試し、局所SPD射影へfallback。エネルギーと0.001mm停止残差は維持。独立非rest応力微分と解析的な均一ひずみ平衡を検証。
- 継続法で未収束試行をaccepted状態へ混入しない。保存NPZとJSONを個別にatomic置換し、snapshotSHAとprogressで混在読みを拒否。旧180°の高コスト計算は既知の境界不良で停止：旧18mm進捗0.475807、12mm0.381976、候補18mm0.603346。完了180°／解像度収束と誤記しない。
- tests/browser-v076-solid.mjsを姿勢と別triangulationに対応。右の静止形状を左にも鏡映する比較器不具合を修正。正式な非対称90°比較はtest-results/v076-conforming-90-side-corrected／v076-conforming-acromial-90-side、各12原画＋2パネル。前者は同じ旧推定境界、後者は起始候補変更を明記。両パネルを目視したが全筋の不良が残り、改善完了ではない。docs/validation/v076/attachment-boundariesへコピー済み。
- Python数値試験16件を全実行し成功。追加は登録・rollback・点cage・PLC・固定自由度細分・global Newton。CIへ追加しrtree1.4.1をoffline requirementsへ追加。最新コミットのActionsは別途確認。Node35はこの継続で再実行し全通過。Web runtimeは未変更。

次は宣言外固定0の最新候補の90°収束と実骨監査を確認し、共通の骨接触拘束を追加する。全角度・全筋、局所体積、自己交差、CPU/GPU・実機、医学的品質は未合格。公開しないことを最終成果とみなさず、実装と検証を続ける。

ユーザーは「引き続き作業して、公開できるレベルまで完成」と再依頼。計画・保存だけで止めず継続する。main／公開は最後の照合でc81256f、今回も変更していない。

- 新規 `tetrahedral_cage.py`／`conforming_surface.py`。元三角形を制御セル境界で分割し、元全頂点と表面を保持する正Jのアフィン変形。三角筋中部14,841/29,678→29,873/59,742、閉鎖辺面数すべて2、元面積相対誤差4.24e−13以下。未分割原型の曲がった三角形を頂点だけで近似する破損を解消できる。
- `constitutive_metric.py`／`projected_solid_metric.py`は厳密応力微分と局所負固有値の0射影。目的関数・固定・J拒否・停止残差は変更しない。独立微分・共有要素・回転の試験を追加。
- `test-results/v076-cage-mid-18-projected/`：coarse cage672点2772要素、conforming source 29,873点。90°外転は全段階収束43.7秒、体積0.9937777、面方向警告0、付着誤差2.85e−14mm。前方90°も全段階収束69.1秒、体積0.9982203、面方向警告0。180°外転は進捗0.44938で上限、残差0.25979mm（未完成）。12mm cageも以前のspline用準備だけで、tetの解像度収束は未実施。
- Actual bone audit `audit-cage-bones.py`：90°外転で肩甲骨に新規深さ増加1.3411mm、上腕骨4.0104mm。全体最大上腕骨深さはrestと同じ4.0753mmだが別頂点が新規貫通。最悪点12473は旧originウェイト0の固定点。
- **新しい重要な発見**：三角筋中部の旧origin固定392点の実肩甲骨面への距離中央値6.683mm、最大19.724mm。最悪貫通点12473のrest gapは10.123mm。322/392点が2mmより遠い。2mmは調査集計のみ、削除基準ではない。旧固定マスクが真の骨付着とは限らず、筋腹／腱膜端と区別が必要。局所固定解除で隠さず、原著モデル・骨面・末端組織から定義データを再設計する。停止側も旧bone-near候補が単なる接触面を含む恐れがある。
- C1 spline試作も作ったが不採用。`v076-spline-mid-18-newton`は90°負荷進捗0.88907でJガードに当たりline-search-stalled、面警告40。geometryは約0.987体積でも局所破綻。spline exact gradient／C1／最小曲げ補完は代数試験を通るが製品合格ではない。
- **注意：GH軸回旋の初回診断はY-X-Yを使っているが、Ludewig2009のGH記述は前方軸の挙上・側方軸の挙上面・上方軸の軸回旋。X-Z-Y等との系列照合がまだ必要。初回診断値をそのままFig5Dへ比較し、骨へ自動補正を加えない。** `audit-axial-sensitivity.py`は探索用・未採用。軸回旋だけで固定点の骨内判定が変わる事実は確認したが、臨床的な基準角の採用は未実装。
- 参考本文再取得：/tmp/v076-ludewig.html と Fig5D /tmp/v076-ludewig-fig5d.jpg、本文Methodsの軸定義と結果・Discussionを確認。現structural scapula frameはTS→AA、論文はTS→AC。図の数値を異なる系列へ直接足さない。
- Seth2019 PMC6856649の本文・筋束表・Data Availabilityを再確認。SimTK project1708 file5812原zipはlogin要求。ログイン不要の同著者らの2023年研究primary repo `ComputationalBiomechanicsLab/rmr-solver` v1.1を参照可能。モデル `OpenSim Models/for RMR solver/TSM_subject_noWeight.osim` とLICENSE_dataを /tmp/v076-reference-model/ に保存。CC BY4.0。BodyPartsへの登録・付着移植は未実施。原論文力推定を本アプリへコピーしない。
- XML DeltoideusScapula_M：origin /scapula [-.0139645,.0019604,.00744838]m、via /scapula[-.00447856,-.00649718,.0296382]、insertion /humerus[-.00524173,-.0851729,.00396206]。このsubject modelのframeであり現world座標ではない。GH CustomJointはplane_elv y、shoulder_elv -x、axial軸[-.0846,.9947,-.0584]。Ludewig系列と同一と仮定しない。

現在の次の一手：固定領域の解剖学的再照合と真の末端／接触の区別、GH座標の照合を先行し、共通接触を加える。conforming cageはソース形状・反転防止の有望な基盤だが、まだ一筋の試作でCPU参照のみ、GPU・インタラクティブ・全筋・公開へ未統合。成功した90°だけで完成扱いしない。
