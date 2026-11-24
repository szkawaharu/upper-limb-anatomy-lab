# v0.61 幾何拘束（研究用・近似）

参考：Seth et al. 2016, A Biomechanical Model of the Scapulothoracic Joint; OpenSim ScapulothoracicJoint https://opensim-org.github.io/opensim-moco-site/docs/1.2.0/html_user/classOpenSim_1_1ScapulothoracicJoint.html ; Seth et al. 2019 https://pmc.ncbi.nlm.nih.gov/articles/PMC6856649/ ; Kavan et al. 2008 https://dcgi.fel.cvut.cz/publications/2008/kavan-tog-adqb . 独自の近似でありOpenSimの移植ではない。

実肋骨表面サンプル範囲から軸平行楕円体を算出し、肩甲骨内側候補点の中立位に対する放射方向距離差を目的関数とする。肩鎖関節中心を固定し2軸各8°以内で正則化付き勾配補正。代理損失が増える補正は採用しない。中立姿勢、同期OFF、胸郭未取得時は従来運動。楕円体は実表面ではなく、真の滑走・貫通防止を保証しない。

筋変形：v0.60.1、v0.60（現固定候補維持）、中間、6回隣接平滑化を候補とし、13姿勢・最大180辺で対数辺長比の二乗と外れ値を評価。最良の固定ウェイトを読み込み時に選択。同じサンプル・姿勢・評価関数では損失が増えないが全姿勢、各辺、筋束長、体積、衝突は保証しない。GPU線形補間を維持、非線形スキニングは未実装。検査UIに代理誤差と候補選択を表示。実メッシュの比較はブラウザ読み込み時に実行し、実機GPU・性能・解剖学精度は未検証。
