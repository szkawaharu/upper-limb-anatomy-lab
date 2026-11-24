# v0.64.0 — 筋束経路エンジン（研究試作）

## 実装
- BodyParts3Dの中立表面を保持。既存付着候補ウェイトを長手パラメータとし、登録横断軸に6本の連続領域ガイドを構築。実測筋束・筋内腱の分割ではない。
- 独立エンジン muscle-fiber-v064.js と V64_MUSCLES 登録表。新筋は骨付着候補、横断軸、包絡種別を登録。未登録・構築不能は従来表示。17解剖筋部位（三角筋3、大胸筋3、前鋸筋、広背筋、僧帽筋3、菱形筋2、肩甲挙筋、肩甲下筋、小胸筋、鎖骨下筋）を登録済み。
- 中立形状、固定ウェイト0/1の剛体変換を維持。局所接線の長手伸縮比λを0.5–2に制限し、横断方向をλ^(-1/2)に。理想局所変形の行列式は1だが、補間・端部テーパ・経路変化のため筋全体の体積保存ではない。
- 三角筋・大胸筋・肩甲下筋は骨頭近傍、前鋸筋・広背筋は胸郭断面の近似包絡。肋骨別の確定筋束や大胸筋腱の捻転は未実装。上腕・前腕の既存模式筋の一部は中立位円錐台体積を参照して端部半径を固定した補正。参照未取得時は旧方式。
- WebGL頂点シェーダとCPU参照で候補を評価。辺長損失、別サンプルの面積損失、該当筋の中立形状相対包絡逸脱が悪化しない候補だけを選択。補間結果も再検査。候補なし・GPU非対応・エンジン取得失敗は従来方式。全頂点の接触・自己交差は検証しない。
- 肩甲骨肋骨面の候補点を優先し、胸郭断面からの中立相対距離を5変数で補正。SC–AC長と接続を維持。実接触・関節包拘束ではない。
- スマホの自動フィット余白を共通化。筋の境界はサンプル点近似であり全頂点包含を保証しない。

## 検証
- node tests/v064-engine.mjs : 合成チューブ7姿勢、中立形状、固定端、有限値、決定性、辺長・面積・包絡代理損失の非悪化、体積代理値。
- 合成肩帯105姿勢（挙上7×挙上面5×軸回旋3）で有限値・鎖骨長・代理損失の非増加を確認。90姿勢で代理損失低下。
- Chromium + SwiftShader WebGL1で頂点/断片シェーダのコンパイル・リンク成功。iOS Safariの描画・性能は未確認。
- tests/v064-atlas-results.json : 固定コミットの実大胸筋腹部・鎖骨部を簡略ランドマークウェイトと上腕のみの回転で比較（本番の肩帯運動・包絡なし）。計10条件で辺長・面積損失非悪化、4条件で新方式採用。三角筋3部位のSTLは取得制限により検証不能。全17部位の実データ検証は未完了。
- window.__labFiberV064() と検査ボタンで左右の現在姿勢・適用率・指標を表示。閉鎖メッシュ体積検査にも対応（自己交差未検証）。

## 文献と採用範囲
1. Seth et al. (2016). A Biomechanical Model of the Scapulothoracic Joint to Accurately Capture Scapular Kinematics during Shoulder Movements. https://doi.org/10.1371/journal.pone.0141028 — 胸郭面上の肩甲骨運動の拘束という考え方。論文モデルの再実装ではない。
2. Ludewig et al. (2009). Motion of the Shoulder Complex During Multiplanar Humeral Elevation. https://pubmed.ncbi.nlm.nih.gov/19181982/ — SC、AC、GHの三次元運動の根拠。個人実測値は未採用。
3. Blemker, Pinsky & Delp (2005). A 3D model of muscle reveals the causes of nonuniform strains in the biceps brachii. https://doi.org/10.1016/j.jbiomech.2004.04.009 — 非一様ひずみ・体積保持を考慮する根拠。有限要素の構成式・活性モデルは未実装。
4. Teran et al. (2005). Creating and Simulating Skeletal Muscle from the Visible Human Data Set. https://graphics.cs.wisc.edu/Papers/2005/TSBNLF05/ — 筋束場、準非圧縮性、有限要素の必要性。体積メッシュは将来課題。
5. OpenSim Muscle Editor. https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim33/pages/53674058/Muscle%2BEditor — 付着、経由点、wrap objectの分離。OpenSimの物理計算そのものではない。
6. Sakoma et al. (2011). Anatomical and functional segments of the deltoid muscle. https://pmc.ncbi.nlm.nih.gov/articles/PMC3042752/ — 三角筋の複数区域。報告の7区域と本版の6本/部位は対応しない。
7. Anatomy, Thorax, Serratus Anterior Muscles. https://www.ncbi.nlm.nih.gov/books/NBK531457/ — 上位8–9肋骨と肩甲骨内側縁への付着。現在の領域ガイドは肋骨別確定筋束ではない。
8. Three-dimensional study of pectoralis major muscle and tendon architecture. https://pubmed.ncbi.nlm.nih.gov/19291757/ — 筋・腱の複雑な構造。筋内腱と停止の詳細は未実装。

## 医学部教育レベルへのリリース条件
- 筋の起始・停止・筋束・筋内腱・骨接触面を解剖学文献と原データで登録し専門家レビュー。
- 左右・全可動域の実データ試験、体積・自己交差・骨/隣接筋貫通・関節面整合の検査。
- 動態MRI/CT・運動学等の独立データと比較。個体差・正常範囲を明示。
- iOS/Android/PC実機性能、入力連続性、アクセシビリティの検証。近似と検証済み部位を区別。
- 形態学習と力/活性の推定を混同しない。医学教育の補助教材として専門家評価を経る。

本版は教育・研究試作。臨床診断・障害リスク・実射負荷・医学部教育適格性を保証しない。
