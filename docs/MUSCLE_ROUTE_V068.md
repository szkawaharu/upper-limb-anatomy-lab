# v0.68.0 — 肩甲帯周囲筋の順序付き経路補正

## 目的
v0.67で肩甲骨・鎖骨・上腕骨の運動が改善した結果、筋メッシュ側に残った以下の問題を減らす。

- 胸郭を直線的にショートカットする筋腹
- 上腕骨頭付近で内側へ巻き込みすぎる三角筋・大胸筋・広背筋
- 肩甲骨運動に追従するときの局所的な折れ・筋腹の急な厚み変化
- 筋束レーン間に生じる段差

## 実装
新規 `muscle-route-v068.js` を追加し、v0.64筋束エンジンが生成した6本の領域ガイド曲線の **中間制御点だけ** を後処理する。起始・停止側の制御点は変更しない。

### 順序付き障害物
OpenSimのgeometry path / wrap objectの考え方を参考に、複数障害物を適用するときは順序を固定した。

1. 胸郭表面proxy
2. 上腕骨頭sphere proxy

腋窩付近では胸郭側の経路を先に整え、その後に局所的な上腕骨頭keepoutをかける。

### 胸郭
BodyParts3Dの実肋骨サンプルからv0.63以降で構築している高さ別楕円断面を再利用。
各ガイド制御点の中立位gapを基準に、内方移動と外方移動を筋別の広い範囲へ制限する。

### 上腕骨頭
v0.62以降の上腕骨表面から得る近似半径を使用。
中立位の半径関係を壊しすぎないよう、筋別の経路区間だけkeepoutを適用。

### 対象
16筋:
- 大・小菱形筋
- 肩甲挙筋
- 肩甲下筋
- 僧帽筋上・中・下部
- 広背筋
- 大胸筋 鎖骨部・胸肋部・腹部
- 小胸筋
- 鎖骨下筋
- 三角筋 前・中・後部

前鋸筋はv0.65の第1–9肋骨9スリップ専用経路を維持し、v0.68の一般経路補正で上書きしない。

## 変形品質
v0.68でガイドを補正した後も、v0.64の既存quality gateを通す。
- 局所辺長損失
- hold-out三角形の面積損失
- 既存の簡略包絡指標

筋束エンジン側では長手方向伸長に対して横断方向を `1/sqrt(lambda)` で補正する近似を維持する。これは局所的な非圧縮近似であり、有限要素法による厳密な体積保存ではない。

## 工学的上限
筋別に、中立位から許す胸郭gap変化、上腕骨頭周囲の外方移動、ガイド制御点の最大移動を設定する。
これらは **臨床正常値や実測筋間距離ではない**。iOS上での破綻防止と、アトラス中立形状から過度に外れないための安全上限。

## 文献・実装上の根拠
1. OpenSim Muscle Editor / Geometry Path
   https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53090145
   - 筋経路にvia point / wrap pointを置き、sphere・ellipsoid・cylinder等の簡略障害物で経路を拘束する考え方。
2. OpenSim FAQ
   https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53084548/Frequently%2BAsked%2BQuestions
   - 任意メッシュではなく簡略wrap objectを使うという実装上の考え方。
3. Lines of action and stabilizing potential of the shoulder musculature. PMCID PMC2740966.
   - 大胸筋・広背筋等のline of actionで、骨性隆起や他筋を回り込むobstacle-set modelが使用されている。
4. Finite Element Modeling of Passive Material Influence on the Deformation and Force Output of Skeletal Muscle. PMCID PMC4075494.
   - 筋をほぼ非圧縮材料として扱うこと、収縮・曲率に伴って局所線維方向が変化すること。
5. A Biomechanical Model of the Scapulothoracic Joint to Accurately Capture Scapular Kinematics during Shoulder Movements. PLoS ONE 2016;11:e0141028.
   - 肩甲骨と胸郭を独立した剛体・曲面拘束として扱うv0.66以降の肩甲帯基盤。

## 起動UXの変更
- `関節操作終了時、画面外の骨格を自動フィット` を初期OFF。
- `姿勢変更時にも自動フィット` は従来どおり初期OFF。
- 実骨読込時のカメラを背面(0 rad)に強制せず、正面寄り斜位から開始。
- 初回だけモデル全体を画面内へ収めるfitは行う。以降はユーザーがフィット操作をONにしない限り、関節操作後に勝手に画角を変えない。

## キャッシュ
`Cache Storage` に選択したBodyParts3D STLを最大48 MiB保存する。
- cache名に固定ソースcommit識別子を含む。
- 1ファイル4 MiB以下のみ対象。
- 読込時に期待byte数を検証。
- 不一致・parse失敗なら該当cacheを削除してネットワーク再取得。
- 保存前にもparse済みデータであることを確認。

したがって表示される「端末キャッシュから再利用」は、古いアプリコードへ戻っていることを意味しない。主目的はモバイル回線での再取得量と待ち時間を減らすこと。

## 検証
`node tests/v068-muscle-route.mjs`
- 16筋 × 5障害物変位 = 80条件
- NaN/Infなし
- 起始・停止制御点を完全固定
- 障害物損失を悪化させない
- 筋別最大ガイド移動を超えない

開発時の合成80条件では57条件で補正が採用された。

## 限界
- 実際の筋同士の接触・押し合いは未実装。
- 肩甲骨・上腕骨・胸郭の任意メッシュに対する完全な衝突判定ではない。
- 筋膜、滑液包、腱鞘、腱滑走、筋力、活動度、pennationの動的変化は未実装。
- 体積保存は近似。
- 医学教育用としては、実機スクリーンショット、独立画像データ、解剖学専門家レビューが必要。
