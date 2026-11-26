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
