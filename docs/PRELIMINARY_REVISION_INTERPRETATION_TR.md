# İlk T1/T2 Sonuçlarının Bilimsel Yorumu

## Mevcut karar

Frozen `structure_transfer_safe + Extra Trees` modeli şimdilik korunmalıdır. İlk sonuçlar modeli
zayıflatmıyor; tersine, CMLase yanıtında gerçek ve siteye doğru eşlenmiş kompakt bir yapısal sinyal
bulunduğunu destekliyor. Ancak bu sinyal yalnız BSA/CrGO-897 deney bağlamı için gösterilmiştir.

## Temel bulgular

1. `structure_full + Extra Trees` MAE=0.1815 ile en iyi sonucu verdi; fold-wise mean null modele göre
   yaklaşık %36.5 iyileşme sağladı.
2. AlphaFold-transfer uyumlu dokuz özellikli `structure_transfer_safe + Extra Trees` MAE=0.1893 ile
   null modele göre yaklaşık %33.8 iyileşme gösterdi.
3. Sequence-only Extra Trees MAE=0.2979 ile null modelden yaklaşık %4.2 daha kötüydü. Sequence
   özellikleri tek başına yararlı değildir.
4. Combined modeller structure-only modellerden belirgin biçimde daha kötüydü. Bu, 101 dizi
   özelliğinin küçük örneklemde gürültü ve varyans eklediğini gösterir.
5. B-factor çıkarılması Extra Trees structure-only modelinde MAE'yi 0.0078 artırdı; katkı olasıdır fakat
   henüz doğrudan paired güven aralığı hesaplanmadığından kesin mekanistik kanıt değildir.
6. Contact proxy combined modelde Extra Trees ve LightGBM için küçük iyileşme, Ridge için kötüleşme
   gösterdi. Sequence gürültüsüyle karıştığı için proxy katkısı henüz izole edilmemiştir.
7. Tek `ward_g4` negatif kontrolünde doğru eşlenmiş yapı MAE=0.2475, permütasyonlu yapı ortalaması
   MAE=0.3157 idi. Doğru eşleme yaklaşık %21.6 daha düşük hata verdi; p=0.004975, 200 iterasyonun
   ulaşabileceği en küçük p-değeridir.

## Henüz söylenmemesi gerekenler

- B-factor'ın aktarılabilir bir esneklik ölçüsü olduğu söylenemez.
- Contact proxy'nin B-factor'ın yerini aldığı söylenemez.
- Proteom çapındaki skorların kalibre edilmiş onarım olasılıkları olduğu söylenemez.
- Tek partisyonlu T2 sonucu bütün spatial bölmelere genellenemez.

## Sonraki karar kapıları

- Proxy, structure-only modelde çoğu partisyonu kazanır ve paired CI sıfırın üstünde kalırsa aktarım
  modeline aday olabilir.
- Nested model seçimi Extra Trees'i çoğu dış fold'da seçer ve sabit Extra Trees nested seçimden kötü
  değilse frozen model korunur.
- T2-full bütün partisyonlarda pozitif delta verirse “structure-aware” iddiası güçlü biçimde korunur.
- Aksi durumda iddia “compact structure-associated signal” biçiminde yumuşatılır.
