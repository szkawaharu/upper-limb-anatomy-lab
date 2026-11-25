# v0.66.0 — 肩甲胸郭3領域面拘束

## 目的
v0.65で残った高挙上時の肩甲骨内側縁・下角の局所的な浮きと、肩甲骨全体の胸郭面への向きの不整合を改善する。見た目を特定姿勢へ固定するのではなく、肩甲骨を剛体として扱い、胸郭面との相対距離・接線整合・鎖骨接続を同時に評価する。

## 実装
- 独立モジュール `scapulothoracic-v066.js` を追加。
- 肩甲骨肋骨面の候補点を高さで **upper / middle / lower** の3領域へ分割。
- 各領域で中立位に対する胸郭表面距離の変化を評価し、RMSだけでなく局所最大離開・最大内方移動も評価。
- lower（内側縁下部〜下角側）は、winging様の局所離開が目立ちやすいため工学的重みを大きく設定。
- 肩甲骨面法線と胸郭断面の楕円近似から得た局所法線の平行性を、**中立位から悪化させない接線整合指標**として追加。
- SC–AC長を変えない。鎖骨は挙上・後退方向の微小補正と長軸回旋を許し、AC位置は常にSCから同一長の鎖骨ベクトルを回転して決定。
- 肩甲骨はACを基準に剛体回転するため、肩甲骨形状そのものは変形しない。
- GH中心は肩甲骨剛体上のランドマークとして維持し、上腕骨の回転原点もそのGH中心を使う。したがって本モデル内では中心位置の離断は生じない。v0.66監査では上腕軸と肩甲骨面の関係を記述値として表示する。
- v0.66ソルバーが取得できない場合はv0.65の肩甲胸郭補正へフォールバック。

## 自由度
v0.66の補正自由度は6:
1. 肩甲骨 内外旋様補正
2. 肩甲骨 前後傾様補正
3. 肩甲骨 上下方回旋様補正
4. 鎖骨 挙上様補正
5. 鎖骨 前後方向補正
6. 鎖骨 長軸回旋補正

これはISBで定義された臨床角度を直接表すものではない。既存アトラス座標上の補正変数である。

## コスト関数
主目的は胸郭との幾何整合。各領域について以下を加算する。
- 中立位相対距離のRMS
- 局所最大離開の超過
- 局所最大内方移動の超過
- 中立位からの胸郭接線整合低下

補正角には小さな正則化を加え、胸郭幾何の改善がない大きな姿勢変更を避ける。最終候補の幾何コストが補正前より悪い場合は採用せず、補正前姿勢を維持する。

**重要:** 6–10 mm等の離開閾値は数値安定化のための工学的設定であり、正常範囲・病的winging診断基準・実際の軟部組織厚を意味しない。

## 文献に基づく設計判断
1. Seth A, Matias R, Veloso AP, Delp SL. A Biomechanical Model of the Scapulothoracic Joint to Accurately Capture Scapular Kinematics during Shoulder Movements. PLoS ONE. 2016;11:e0141028. https://doi.org/10.1371/journal.pone.0141028
   - 胸郭を曲面として扱い、肩甲骨をその表面上で拘束する剛体モデルの考え方を採用。本実装は同モデルの数式・パラメータの複製ではない。
2. Scapular kinematics variability in individuals with and without rotator cuff-related shoulder pain: systematic review with multilevel meta-regression. PMID 40845626.
   - 健常者でも肩甲骨姿勢・運動には幅があり、挙上に伴い上方回旋・後傾が増え、外旋方向の傾向がある。個人差が大きいため、v0.66では単一の「正常角度」への強制を避ける。
3. Spartacus: A review and aggregation of reference datasets reporting the normal shoulder girdle kinematics during uniplanar humerus motions. PMID 40554241.
   - 研究間で座標系・算出法に不整合があることを踏まえ、臨床角度の固定値を直接コストに使用しない。
4. Mechanistic insights into glenohumeral kinematics derived from positional relationship between the contact path and humeral tuberosity. PMID 36701958.
   - GH運動は単純な球関節固定回転だけではなく接触経路を伴う。本版はGH中心を固定した教育用剛体近似であり、接触経路・軟骨・回旋腱板による求心化を再現しない。
5. Determination of inter-individual differences in humeral head centre translations in subjects with healthy shoulders. PMID 42841165.
   - 健常肩でもHHC translationには個体差と運動依存性がある。v0.66のGH中心固定は計算上の近似であり、生体内translationの再現ではない。

## 検証
- `node tests/v066-scapulothoracic.mjs`
  - 7挙上角 × 5挙上面 × 3軸回旋条件 = 105条件
  - 有限値
  - SC–AC長保持
  - 3領域サンプル成立
  - 幾何目的関数の非悪化
- JavaScript構文検査。
- アプリ内 `v0.66 肩甲胸郭・GH監査` から、upper / middle / lower のRMS・最大離開・最大内方・接線整合を現在姿勢で表示可能。
- `window.__labScapulothoracicV066()` で同監査値を取得可能。

## 未検証・今後必要
- iPad/iPhone実機で全方向・全可動域の見た目と性能。
- 動態CT / biplanar fluoroscopy / bone-pin等の独立データとの定量比較。
- 肩甲骨と胸郭の実接触面、前鋸筋・肩甲下筋・滑液包を含む軟部組織厚。
- AC/SC/GHの関節面・靱帯・関節包・筋力による拘束。
- GH contact pathと上腕骨頭translation。
- 左右差・年齢・性別・個体差を扱うパラメータ化。
- 解剖学専門家によるレビュー。

本版は教育・研究試作であり、臨床診断、正常/異常判定、障害リスク評価には使用しない。
