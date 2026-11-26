# 未採用の筋変形実験

これらをindex.htmlからロードしない。テスト通過は合成入力に対する代数・境界・状態分離の検証であり、人体形状や医学的正確性の合格ではない。[採否](../../docs/MUSCLE_DEFORMATION_V076.md)と[最新引継ぎ](../../WORKING_STATE_V076.md)を読む。

## データを再生成する順序

固定アトラスコミットは2c83c698788cff57f6317bddcf712e7acfbea7d8。STLキャッシュは/tmp/upper-limb-atlas。Python依存はscripts/calibration/requirements.txtを参照（実行環境のPython例：/tmp/upper-limb-research-venv/bin/python）。コマンドのpythonはその環境に置き換える。

```sh
QA_NUMERIC_ONLY=1 QA_EXPORT_BINDINGS=1 QA_SOURCE_MESH=1 QA_OUTPUT=test-results/v076-source-bindings node tests/browser-v073.mjs
python scripts/calibration/restore-source-topology.py
python scripts/calibration/audit-rest-attachments.py
python scripts/calibration/infer-terminal-regions.py
python scripts/calibration/build-harmonic-attachments.py test-results/v076-source-terminals/bindings.json test-results/v076-source-attached
node scripts/calibration/refit-source-guides.mjs
python scripts/calibration/build-muscle-thickness.py test-results/v076-source-clean/bindings.json test-results/v076-source-clean/thickness
python scripts/calibration/build-source-material.py test-results/v076-source-refit/bindings.json test-results/v076-source-refit/material
python scripts/calibration/build-deformation-graphs.py
python scripts/calibration/build-source-volume-quadrature.py
```

`restore-source-topology.py`はexact weldのfirst-appearance番号と元STLが一致することをassertする。違うアトラスやクラスタ済みexportへ強行適用しない。source/exportのSHAを保存し、入力が変わった場合は別出力ディレクトリにする。

`infer-terminal-regions.py`は未取得の固定端だけをsource骨面と照合する。ROIも候補であって臨床確認済みの付着ではない。`build-deformation-graphs.py`には、大胸筋腹部の未照合起始に対する**低信頼のinferior terminal-cap仮説**が一箇所ある。原型の下端を体幹へ固定する実験用境界条件であり、腱膜の実測付着と呼ばず製品へ持ち込まない。新規筋の共通デフォルトにしてはならない。

元のクラスタ済み57e3530のbindingはtest-results/v076-material-baseline。再生成した新入力で以前の筋名・頂点番号・数値の完全再現を主張しない。exportは実行したコードに依存するため、毎回source-manifest.jsonを残す。

## ARAP参照と反例

```sh
QA_GRAPH_VOLUME=1 QA_GRAPH_VOLUME_DIR=test-results/v076-source-graphs/volume QA_GRAPH_CONTACT=1 QA_GRAPH_ELEVATIONS=0,90,180 QA_GRAPH_TAG=exact-mid-final QA_GRAPH_FILTER=v60DeltMid node scripts/experimental/audit-deformation-graph.mjs
node tests/v076-deformation-graph.mjs
```

- `QA_GRAPH_VOLUME_DIR`未指定時は旧層別求積を使用。求積近似が正確な全source体積から外れる場合があり、両者を区別する。
- `QA_GRAPH_BOUNDARY=1`は未採用の四次境界基底、`QA_GRAPH_INITIALIZER=dqs`は初期値比較、`QA_GRAPH_EMBEDDING=orthogonal`は回転平均比較。比較条件を出力manifestへ記録。
- 接触ONは校正済み骨頭球の節点近似。全骨の厳密な三角形衝突ではない。
- 全体積の制約が解けても、面方向警告・自己交差・骨貫通・性能の合格を意味しない。

## 他の未採用案

`rest-surface*.js`は微分面積と半弦による補正、`rest-transport.js`は中心／直線／DQSとスカラー体積投影の比較。半弦レイは筋厚さの実測ではない。`surface-transport-volume-v1.js`等は途中スナップショット。`audit-chart-candidates.mjs`の入力名は途中段階を指すため、そのファイルが存在する条件でのみ再実行する。

`build-shoulder-field-uniform.mjs`は不採用の均一格子実験。採用済みの肩甲帯補間を生成する場合はscripts/calibration/build-shoulder-field.mjsを使い、1mmの判定を維持する。

`pack-source-bindings.mjs`とtests/helpers/source-preview.mjsは原型をCPU描画するテスト専用。QA_PREVIEW_SOURCE=1で明示した比較だけに使用し、製品版やGPU変形の比較と混同しない。

## 原型に一致する四面体参照（試作）

[参照計算の状態](../../docs/VOLUMETRIC_REFERENCE_V076.md)。既定へ統合していない。Python研究環境でtetgen==0.8.4を別途使用する。optional生成ツールbuild-volumetric-reference.pyはAGPL-3.0-or-later、同梱したcoreのライセンスはlicenses/tetgen-license.txt。Webアプリへこの生成ツールやTetGenバイナリを同梱せず、既存ファイルのライセンスも変更しない。元の解剖データの出典は保持する。

```sh
python scripts/experimental/build-volumetric-reference.py
node scripts/experimental/export-volumetric-targets.mjs
python scripts/experimental/test-volumetric-reference.py
python scripts/experimental/run-volumetric-reference.py test-results/v076-volumetric 0 0 100
```

生成は境界の修復・簡略化をしない。元の全表面頂点・面の完全一致と元の閉鎖体積との一致を検査し、失敗したら停止する。最初のdriverは監査済み三角筋中部に限定。参照データを再生成する場合は過去の解を上書きせず別ディレクトリへ出す。局所の体積外れ値・停止量・骨接触・性能の問題が残り、非反転だけで完成としない。

## 体積参照の追加解法と軽量化の反例

```sh
python scripts/experimental/test-volumetric-gn.py
python scripts/experimental/test-volumetric-transfer.py
python scripts/experimental/test-reduced-volume.py
python scripts/experimental/test-elastic-modes.py
python scripts/experimental/refine-volumetric-gn.py INPUT.npz 3000 30 OUTPUT
# 連成source-rest前処理（変形の許容値や物性は変えない）
python scripts/experimental/refine-volumetric-gn.py INPUT.npz 3000 30 OUTPUT coupled
python scripts/experimental/prepare-reduced-volume.py test-results/v076-volumetric-quality 24
python scripts/experimental/prepare-elastic-modes.py test-results/v076-volumetric-quality 30
python scripts/experimental/run-vector-continuation.py test-results/v076-volumetric-quality elastic-modes-k30 0 0 30 40
```

`requirements-reference.txt`はTetGenなしの代数テスト用。体積分割の生成環境とは区別する。`build-volumetric-reference.py`のquality／出力引数は実コードを確認して別ディレクトリへ出す。

`reduced_volume_basis.py`のsource 4チャンネルFEM調和liftingは、両方の元固定領域と共通剛体運動を厳密に保つ。スカラー場を筋活動と解釈しない。`linear_elastic_modes.py`は元のエネルギーの厳密なrest接線応答を使うが、`vector_reduced_solid.py`の少ない自由度が全FEMを再現したとは限らない。結果には射影残差と全FEM残差を別々に出す。どちらの縮約案も実メッシュ中立への継続で局所破綻し、不採用。途中停止したNPZを中立位の完成結果として表示しない。

現在変形での連成前処理は `refine-volumetric-gn.py INPUT.npz 3000 100 OUTPUT frozen`。不採用の `coupled`（restだけ）と `warped`（節点回転だけ）から区別する。エネルギー・固定端・停止閾値は同一。独立監査は `audit-solid-reference.py INPUT.npz OUTPUT.json`。

```sh
python scripts/experimental/prepare-solid-preview.py test-results/v076-source-refit/bindings.json test-results/v076-volumetric-quality/neutral-refined-gn-k3000-frozen-v2.npz test-results/v076-volumetric-quality/neutral-refined-gn-k3000-frozen-v2.json test-results/v076-volumetric-quality/neutral-preview.json.gz
QA_SOLID_INPUT=test-results/v076-volumetric-quality/neutral-preview.json.gz QA_OUTPUT=test-results/v076-solid-quality-neutral node tests/browser-v076-solid.mjs
python scripts/experimental/build-solid-panels.py test-results/v076-solid-quality-neutral
```

previewの生成が完了してからブラウザを起動する。元の全表面点・面・固定領域と、宣言した中立の骨境界が合わなければexportを停止する。静止snapshotをインタラクティブな変形の実装完了と呼ばない。パネル生成のみ追加でmatplotlibが必要（研究venv）、CIの5代数試験には不要。
