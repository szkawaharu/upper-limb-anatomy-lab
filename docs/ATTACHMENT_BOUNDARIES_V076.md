# v0.76：付着領域と補間固定範囲の修正（未採用）

2026-11-27。製品runtimeは変更していない。公開はv0.75.0のまま。骨格の最終姿勢を変えず、元の解剖表面と固定境界の意味を再監査している。

## 確認した原因

1. 旧「起始固定」の推定に骨付着以外の組織が含まれる。中部三角筋392点の肩甲骨面への距離は中央値6.683mm、最大19.724mm。90°外転の新規上腕骨内深さ4.0104mmは旧固定点12473自身で、筋変形だけでは解けない。
2. 粗い格子では付着点の補間支持点をすべて固定するため、宣言していない周囲の組織まで固定される。新しい起始候補でも4,313表示頂点に自由な支持点がなく、最悪点28359はラベルが約0.098でも完全に固定されていた。
3. 付着頂点だけをDelaunay制御点へ追加しても、付着三角形の内部が自由制御点へ補間される場合がある。合成の点拘束試験だけでは発見できず、実データの面拘束検査で拒否した。
4. PLCで付着面を回復しても、実際の付着ではない内部の制御辺／面が、固定点だけで構成されると剛体化する。残る57点を安静時の表面と補間から特定した。

## 実装した構造

- 起始候補の定義を `scripts/experimental/acromial-region-candidate.json` へ分離。肩峰外側壁の解剖区分と元骨メッシュから抽出し、起始は497点（旧392点との共通57点）。元の頂点・面は変更しない。骨面距離0.00049〜1.54715mmは候補抽出の結果で、医学的合格値ではない。帯幅1.55108mmは測定した骨／筋メッシュの辺長から算出し、姿勢や見た目を評価して調整していない。腱膜・筋膜を含む実付着範囲は未確認。
- PLC制御メッシュで、起始／停止の元993頂点と1,702固定面を保持。通常の点Delaunayを採用しない。TetGenはAGPLのオフライン生成器に限定し、Webへリンクしない。
- `attachment_refinement.py` と `refine-attachment-supports.py` は、宣言していない表面組織の剛体化を安静時だけで検出し、該当制御辺・面・セルへ自由な制御点を追加する。17辺、80面、84セルに計181点を追加。元の付着面・制御点を動かさず、親セルの積分密度を体積比で保存する。
- 制御点2,863→3,044、要素17,410→18,124。表示38,712点の未宣言剛体化0。元14,841点と29,678三角形の表面を保持し、表示77,420面へ分割。各元三角形の面積和の誤差最大4.57e−13、全116,130辺の接続面数2。
- 数値負荷の未収束試行をロールバックする継続法を追加。目標姿勢・エネルギー・固定端・0.001mm停止残差は変更しない。保存は一時ファイルからの置換とSHA照合で、途中のZIPや異なる段階の組を拒否する。
- 厳密な全体接線のNewton方向を、下降方向と正J・Armijo条件を満たす場合に利用。局所固有値射影の方式へフォールバックできる。材料定数や許容値を変えた高速化ではない。

## 比較と未合格項目

[検証データ・比較画像](validation/v076/attachment-boundaries/)を保存。90°外転の曲線方式／制御格子方式を、同じ元表面・最終骨・カメラで正面／側面／背面、全筋／右中部三角筋の計12画像ずつ取得した。起始候補変更あり／なしを別記した。ブラウザ比較が右の静止形状を左へも鏡映する不具合を修正し、左右非対称の画像を撮り直した。初回の誤った非対称画像は正式資料へコピーしない。

未細分PLCの90°外転は全段階収束、体積比0.999074、固定点誤差0、面方向警告0。ただし新規上腕骨内深さ2.27577mm、辺長比0.303〜3.199、面積比0.132〜3.269で合格ではない。起始候補＋粗格子の4.00902mmより浅くなったことだけを採用理由にしない。新しい未宣言固定0のPLCは解き直し中。接触拘束、全可動域、全筋、局所体積、自己交差、GPU・実機性能は未達。

旧マスクの180°負荷は18mm格子で進捗0.475807、12mmで0.381976、起始候補＋18mm格子で0.603346までの途中結果。既知の境界条件不良で計算を止めたもので、180°の合格や解像度収束比較とは呼ばない。

外部モデルは[著者らの公開v1.1データ](https://github.com/ComputationalBiomechanicsLab/rmr-solver/tree/v1.1/OpenSim%20Models)を取得・登録試行したが、TS/AA/AIの剛体相似対応に7.8〜12.3mmの残差がある。モデルAAマーカーも表示骨面から15.38mm離れており、実骨ランドマークと同一とは扱えない。座標・筋経路の移植は未採用。データのCC BY 4.0と著者表記を保持する。

肩峰外側壁の広い腱性付着と上方の筋膜への連続性は、[Smithら2020、本文・組織所見](https://doi.org/10.1016/j.asmr.2020.06.014)を確認した。三角筋を単一の直線とみなさない根拠として、[Sakomaら、60肩の解剖研究](https://doi.org/10.1111/j.1469-7580.2010.01325.x)の筋内腱と区分を確認。これらの群データを当アトラス個体の付着実測や材料定数として扱っていない。

## 再生成

Pythonは研究venv（numpy/scipy/rtree、生成と実骨監査には別途TetGen/trimesh/embreex）を使う。下記は `python` をその環境へ置換する。入力bindingsと全原型四面体の生成は既存の実験READMEを参照。

```sh
python scripts/experimental/prepare-acromial-candidate.py test-results/v076-source-refit/bindings.json test-results/v076-acromial-candidate
python scripts/calibration/build-harmonic-attachments.py test-results/v076-acromial-candidate/bindings.json test-results/v076-acromial-harmonic
python scripts/experimental/prepare-boundary-volume-input.py test-results/v076-volumetric test-results/v076-acromial-harmonic/bindings.json test-results/v076-acromial-volume-input
python scripts/experimental/prepare-constrained-cage.py test-results/v076-acromial-volume-input 18 test-results/v076-plc-mid-acromial-18
python scripts/experimental/refine-attachment-supports.py test-results/v076-plc-mid-acromial-18 test-results/v076-plc-surface-refined-mid-acromial-18
python scripts/experimental/audit-fixed-support.py test-results/v076-plc-surface-refined-mid-acromial-18
python scripts/experimental/run-conforming-reference.py test-results/v076-plc-surface-refined-mid-acromial-18 90 0 30
python scripts/experimental/audit-cage-bones.py test-results/v076-plc-surface-refined-mid-acromial-18/reference-e90-p0-k30.npz
```

元17,410要素のPLCから細分後へ積分密度を継承する。`prepare-constrained-cage.py` の `full-refinement` は全域細分（11,332点／54,825要素）の比較用で、軽量版の標準生成手順ではない。
