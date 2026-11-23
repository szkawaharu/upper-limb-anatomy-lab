# 上肢3D解剖ラボ（Upper Limb 3D Anatomy Lab）

スマートフォン対応の上肢骨格・筋肉・腱を観察する教育用Webアプリ。動作と形態は**解剖学的に未検証の近似を含み、臨床用・診断用ではありません**。

**現在のバージョン：v0.53**

- **アプリ**：https://szkawaharu.github.io/upper-limb-anatomy-lab/
- **公開ページのコード**：[`index.html`](./index.html)
- **今回の肩甲骨筋修正／制限**：[`docs/SCAPULAR_MUSCLES_V053.md`](./docs/SCAPULAR_MUSCLES_V053.md)
- **前版の記録**：[`docs/SCAPULAR_MUSCLES_V051.md`](./docs/SCAPULAR_MUSCLES_V051.md)
- **v0.51の保存ブランチ**：[`archive/v0.51-before-scapular-fix`](../../tree/archive/v0.51-before-scapular-fix)

## 主な機能

BodyParts3Dの肩甲骨、鎖骨、上腕骨、橈骨、尺骨、および右手の骨27個を利用した表示（データ取得不可時は近似形状）。肘屈伸・肩挙上・前腕回内外・手指の教育用の概形操作を提供します。筋肉・腱は手作業の推定パスです。

v0.53では棘上筋・棘下筋・小円筋・大円筋の**後面投影と筋腹の厚み**を修正。小円筋までの3筋は回旋筋腱板、大円筋は含まれません。**肩甲下筋は未実装**。

## 開発・テスト

HTML/JSのみで起動しますが、実骨メッシュを表示するにはオンライン通信が必要です。

```bash
node tests/scapular-muscles-v053.mjs
node tests/scapular-muscles-v051.mjs
node tests/medical-anatomy-audit.mjs
node tests/shoulder-kinematics-v050.mjs
```

テスト通過は3Dの骨貫通・生体組織の滑走・医学的な位置精度を保証しません。実骨データ上の筋の起始・停止・腱経路の照合と専門家レビューが必要です。
