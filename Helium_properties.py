import numpy as np
from scipy.interpolate import interp2d
#monophase helium properties tables obtained from HEPAK (liquid or gas)

P = np.array([0.014,0.016,0.02,0.05,0.1,0.4,0.8,1.0,1.2,1.25,1.3,1.4,1.5,2.0,2.5,3.0,5.0,7.0,10.0,12.0,15.0,17.0,19.0])
T = np.array([1.60,1.70,1.80,2.00,3.00,4.00,5.00,6.00,7.00,8.00,9.00,10.00,20.00,30.00,40.00,50.00,60.00,70.00,80.00,90.00,100.00,110.00,120.00,130.00,140.00,150.00,160.00,170.00,180.00,190.00,200.00,210.00,220.00,230.00,240.00,250.00,260.00,270.00,280.00,290.00,300.00,310.00,320.00,330.00,340.00,350.00])

#density
D_PT = np.array([[145.29,145.29,145.30,145.36,145.45,145.98,146.67,147.01,147.35,147.43,147.51,147.67,147.84,148.64,149.42,150.19,153.07,155.72,159.33,161.53,164.58,166.49,168.32],
[145.34,145.35,145.35,145.41,145.50,146.04,146.74,147.08,147.42,147.51,147.59,147.76,147.92,148.73,149.53,150.30,153.22,155.90,159.55,161.78,164.88,166.83,168.70],
[0.39,0.44,145.42,145.48,145.57,146.12,146.83,147.18,147.53,147.61,147.70,147.87,148.03,148.86,149.67,150.46,153.42,156.14,159.84,162.10,165.29,167.30,169.26],
[0.34,0.40,0.50,145.69,145.79,146.37,147.12,147.49,147.85,147.94,148.03,148.21,148.39,149.27,150.12,150.94,154.03,156.86,160.77,163.24,166.84,168.96,170.86],
[0.23,0.26,0.32,0.82,1.70,141.79,142.83,143.33,143.82,143.94,144.06,144.29,144.53,145.66,146.73,147.74,151.41,154.58,158.69,161.13,164.43,166.45,168.35],
[0.17,0.19,0.24,0.61,1.24,5.45,13.16,129.86,130.81,131.04,131.26,131.70,132.13,134.12,135.90,137.51,142.89,147.17,152.42,155.40,159.35,161.71,163.90],
[0.14,0.15,0.19,0.49,0.98,4.13,8.97,11.78,14.97,15.84,16.75,18.70,20.86,102.85,112.39,117.53,128.98,135.89,143.29,147.20,152.16,155.03,157.65],
[0.11,0.13,0.16,0.40,0.81,3.35,7.04,9.03,11.14,11.69,12.25,13.39,14.57,21.26,29.88,42.21,102.08,118.69,130.89,136.43,142.97,146.59,149.80],
[0.10,0.11,0.14,0.35,0.69,2.83,5.84,7.43,9.06,9.48,9.90,10.76,11.63,16.26,21.43,27.27,59.44,90.96,114.11,122.65,131.74,136.45,140.49],
[0.08,0.10,0.12,0.30,0.60,2.46,5.02,6.34,7.70,8.04,8.39,9.08,9.79,13.45,17.34,21.51,41.34,64.62,93.31,105.80,118.48,124.69,129.84],
[0.07,0.09,0.11,0.27,0.54,2.17,4.41,5.55,6.72,7.01,7.31,7.90,8.50,11.56,14.75,18.08,32.87,49.72,74.62,88.36,103.93,111.73,118.12],
[0.07,0.08,0.10,0.24,0.48,1.95,3.94,4.95,5.97,6.23,6.49,7.01,7.53,10.18,12.91,15.72,27.76,40.95,61.35,73.94,89.98,98.69,106.06],
[0.03,0.04,0.05,0.12,0.24,0.96,1.93,2.41,2.89,3.01,3.13,3.37,3.61,4.82,6.02,7.23,12.05,16.85,24.02,28.75,35.75,40.33,44.82],
[0.02,0.03,0.03,0.08,0.16,0.64,1.28,1.60,1.92,2.00,2.08,2.24,2.40,3.19,3.99,4.78,7.93,11.05,15.66,18.69,23.18,26.12,29.03],
[0.02,0.02,0.02,0.06,0.12,0.48,0.96,1.20,1.44,1.50,1.56,1.68,1.80,2.39,2.99,3.58,5.94,8.27,11.72,14.00,17.36,19.57,21.75],
[0.01,0.02,0.02,0.05,0.10,0.38,0.77,0.96,1.15,1.20,1.25,1.34,1.44,1.92,2.39,2.87,4.76,6.63,9.40,11.22,13.93,15.71,17.47],
[0.01,0.01,0.02,0.04,0.08,0.32,0.64,0.80,0.96,1.00,1.04,1.12,1.20,1.60,1.99,2.39,3.97,5.53,7.85,9.38,11.65,13.14,14.63],
[0.01,0.01,0.01,0.03,0.07,0.27,0.55,0.69,0.82,0.86,0.89,0.96,1.03,1.37,1.71,2.05,3.40,4.75,6.74,8.06,10.02,11.31,12.59],
[0.01,0.01,0.01,0.03,0.06,0.24,0.48,0.60,0.72,0.75,0.78,0.84,0.90,1.20,1.50,1.80,2.98,4.16,5.91,7.07,8.79,9.93,11.06],
[0.01,0.01,0.01,0.03,0.05,0.21,0.43,0.53,0.64,0.67,0.69,0.75,0.80,1.07,1.33,1.60,2.65,3.70,5.26,6.30,7.83,8.85,9.86],
[0.01,0.01,0.01,0.02,0.05,0.19,0.38,0.48,0.58,0.60,0.62,0.67,0.72,0.96,1.20,1.44,2.39,3.34,4.74,5.68,7.07,7.99,8.90],
[0.01,0.01,0.01,0.02,0.04,0.17,0.35,0.44,0.52,0.55,0.57,0.61,0.66,0.87,1.09,1.31,2.17,3.04,4.32,5.17,6.44,7.27,8.11],
[0.01,0.01,0.01,0.02,0.04,0.16,0.32,0.40,0.48,0.50,0.52,0.56,0.60,0.80,1.00,1.20,1.99,2.78,3.96,4.74,5.91,6.68,7.45],
[0.01,0.01,0.01,0.02,0.04,0.15,0.30,0.37,0.44,0.46,0.48,0.52,0.55,0.74,0.92,1.11,1.84,2.57,3.66,4.38,5.46,6.18,6.89],
[0.00,0.01,0.01,0.02,0.03,0.14,0.27,0.34,0.41,0.43,0.45,0.48,0.51,0.69,0.86,1.03,1.71,2.39,3.40,4.08,5.08,5.74,6.41],
[0.00,0.01,0.01,0.02,0.03,0.13,0.26,0.32,0.38,0.40,0.42,0.45,0.48,0.64,0.80,0.96,1.60,2.23,3.18,3.81,4.74,5.37,5.99],
[0.00,0.00,0.01,0.02,0.03,0.12,0.24,0.30,0.36,0.38,0.39,0.42,0.45,0.60,0.75,0.90,1.50,2.09,2.98,3.57,4.45,5.04,5.62],
[0.00,0.00,0.01,0.01,0.03,0.11,0.23,0.28,0.34,0.35,0.37,0.40,0.42,0.57,0.71,0.85,1.41,1.97,2.81,3.36,4.19,4.75,5.29],
[0.00,0.00,0.01,0.01,0.03,0.11,0.21,0.27,0.32,0.33,0.35,0.37,0.40,0.53,0.67,0.80,1.33,1.86,2.65,3.18,3.96,4.49,5.00],
[0.00,0.00,0.01,0.01,0.03,0.10,0.20,0.25,0.30,0.32,0.33,0.35,0.38,0.51,0.63,0.76,1.26,1.76,2.51,3.01,3.76,4.25,4.75],
[0.00,0.00,0.00,0.01,0.02,0.10,0.19,0.24,0.29,0.30,0.31,0.34,0.36,0.48,0.60,0.72,1.20,1.68,2.39,2.86,3.57,4.04,4.51],
[0.00,0.00,0.00,0.01,0.02,0.09,0.18,0.23,0.27,0.29,0.30,0.32,0.34,0.46,0.57,0.69,1.14,1.60,2.28,2.73,3.40,3.85,4.30],
[0.00,0.00,0.00,0.01,0.02,0.09,0.17,0.22,0.26,0.27,0.28,0.31,0.33,0.44,0.55,0.66,1.09,1.52,2.17,2.61,3.25,3.68,4.11],
[0.00,0.00,0.00,0.01,0.02,0.08,0.17,0.21,0.25,0.26,0.27,0.29,0.31,0.42,0.52,0.63,1.04,1.46,2.08,2.49,3.11,3.52,3.93],
[0.00,0.00,0.00,0.01,0.02,0.08,0.16,0.20,0.24,0.25,0.26,0.28,0.30,0.40,0.50,0.60,1.00,1.40,1.99,2.39,2.98,3.38,3.77],
[0.00,0.00,0.00,0.01,0.02,0.08,0.15,0.19,0.23,0.24,0.25,0.27,0.29,0.38,0.48,0.58,0.96,1.34,1.91,2.30,2.86,3.24,3.62],
[0.00,0.00,0.00,0.01,0.02,0.07,0.15,0.19,0.22,0.23,0.24,0.26,0.28,0.37,0.46,0.55,0.92,1.29,1.84,2.21,2.76,3.12,3.48],
[0.00,0.00,0.00,0.01,0.02,0.07,0.14,0.18,0.21,0.22,0.23,0.25,0.27,0.36,0.45,0.53,0.89,1.24,1.77,2.13,2.65,3.00,3.35],
[0.00,0.00,0.00,0.01,0.02,0.07,0.14,0.17,0.21,0.21,0.22,0.24,0.26,0.34,0.43,0.52,0.86,1.20,1.71,2.05,2.56,2.90,3.24],
[0.00,0.00,0.00,0.01,0.02,0.07,0.13,0.17,0.20,0.21,0.22,0.23,0.25,0.33,0.41,0.50,0.83,1.16,1.65,1.98,2.47,2.80,3.13],
[0.00,0.00,0.00,0.01,0.02,0.06,0.13,0.16,0.19,0.20,0.21,0.22,0.24,0.32,0.40,0.48,0.80,1.12,1.60,1.92,2.39,2.71,3.02],
[0.00,0.00,0.00,0.01,0.02,0.06,0.12,0.16,0.19,0.19,0.20,0.22,0.23,0.31,0.39,0.47,0.77,1.08,1.55,1.85,2.31,2.62,2.93],
[0.00,0.00,0.00,0.01,0.02,0.06,0.12,0.15,0.18,0.19,0.20,0.21,0.23,0.30,0.38,0.45,0.75,1.05,1.50,1.80,2.24,2.54,2.84],
[0.00,0.00,0.00,0.01,0.01,0.06,0.12,0.15,0.17,0.18,0.19,0.20,0.22,0.29,0.36,0.44,0.73,1.02,1.45,1.74,2.17,2.46,2.75],
[0.00,0.00,0.00,0.01,0.01,0.06,0.11,0.14,0.17,0.18,0.18,0.20,0.21,0.28,0.35,0.42,0.71,0.99,1.41,1.69,2.11,2.39,2.67],
[0.00,0.00,0.00,0.01,0.01,0.06,0.11,0.14,0.16,0.17,0.18,0.19,0.21,0.27,0.34,0.41,0.69,0.96,1.37,1.64,2.05,2.32,2.59]])

#enthalpy
h_PT = np.array([[396.79,398.18,400.94,421.69,456.26,663.30,938.35,1075.47,1212.32,1246.49,1280.64,1348.90,1417.09,1757.12,2095.61,2432.65,3767.23,5082.03,7021.02,8293.63,10175.98,11415.30,12643.93],
[585.44,586.82,589.59,610.37,644.99,852.37,1127.95,1265.38,1402.56,1436.82,1471.06,1539.50,1607.89,1948.97,2288.68,2627.04,3967.63,5288.71,7237.32,8517.21,10413.39,11664.68,12908.02],
[24262.81,24210.78,844.66,865.49,900.18,1108.02,1384.34,1522.19,1659.83,1694.21,1728.57,1797.27,1865.92,2208.43,2549.74,2889.84,4238.20,5567.60,7530.10,8821.36,10739.93,12010.83,13278.36],
[25386.53,25347.08,25267.80,1655.33,1690.34,1900.19,2179.55,2319.06,2458.45,2493.28,2528.10,2597.72,2667.31,3014.79,3361.44,3707.22,5081.23,6442.28,8470.51,9823.71,11886.86,13147.51,14336.49],
[30736.61,30720.03,30686.79,30434.72,30002.77,5284.47,5509.56,5622.66,5736.07,5764.46,5792.88,5849.75,5906.69,6192.19,6478.76,6766.12,7919.51,9073.57,10797.71,11940.06,13641.22,14766.72,15885.29],
[35971.49,35960.93,35939.79,35780.28,35510.54,33764.84,30830.98,8952.01,9014.57,9030.96,9047.63,9081.75,9116.84,9304.64,9508.56,9724.34,10660.70,11660.94,13214.28,14266.45,15854.60,16915.41,17975.68],
[41183.25,41175.38,41159.63,41041.08,40841.82,39596.67,37762.48,36737.79,35606.95,35302.65,34987.81,34321.21,33593.46,15682.32,14816.82,14554.92,14657.22,15280.76,16505.80,17408.10,18823.87,19791.88,20771.40],
[46387.28,46380.99,46368.39,46273.69,46115.06,45140.32,43768.71,43045.13,42291.03,42097.16,41900.97,41501.25,41091.03,38844.24,36121.42,32586.82,21422.07,20360.83,20750.32,21359.53,22481.24,23305.71,24167.14],
[51587.79,51582.56,51572.10,51493.52,51362.15,50562.35,49461.77,48894.83,48315.45,48168.53,48020.75,47722.50,47420.56,45849.51,44159.33,42326.57,33930.14,28452.02,26411.78,26383.18,26962.00,27550.63,28230.29],
[56786.35,56781.91,56773.02,56706.30,56594.87,55920.33,55003.67,54537.29,54065.13,53946.15,53826.78,53586.88,53345.37,52112.46,50833.97,49506.43,43776.93,38444.18,33926.89,32792.85,32435.58,32643.53,33044.22],
[61983.68,61979.85,61972.18,61914.65,61818.64,61239.49,60458.48,60063.97,59666.67,59566.90,59466.94,59266.48,59065.28,58047.79,57010.57,55953.31,51562.68,47247.47,42327.67,40367.08,38931.80,38630.44,38648.89],
[67180.17,67176.83,67170.14,67119.99,67036.34,66532.87,65857.27,65517.56,65176.54,65091.07,65005.53,64834.18,64662.50,63798.97,62926.95,62046.77,58468.98,54942.87,50419.11,48218.43,46108.88,45344.14,44959.07],
[119126.59,119125.51,119123.34,119107.08,119079.99,118917.87,118702.84,118595.82,118489.13,118462.51,118435.91,118382.77,118329.72,118065.74,117803.93,117544.33,116529.26,115554.52,114178.06,113324.61,112152.38,111448.63,110810.54],
[171062.49,171062.16,171061.50,171056.54,171048.28,170998.93,170933.66,170901.25,170869.00,170860.96,170852.93,170836.90,170820.91,170741.53,170663.14,170585.74,170286.29,170003.65,169612.81,169375.37,169055.58,168867.52,168700.21],
[222995.63,222995.66,222995.72,222996.16,222996.89,223001.39,223007.62,223010.84,223014.12,223014.95,223015.78,223017.46,223019.16,223027.91,223037.09,223046.70,223089.57,223139.79,223229.53,223299.40,223419.90,223511.04,223611.07],
[274927.67,274927.90,274928.37,274931.85,274937.65,274972.50,275019.07,275042.39,275065.74,275071.58,275077.43,275089.12,275100.82,275159.44,275218.23,275277.22,275515.07,275756.15,276124.22,276374.14,276756.21,277015.93,277279.80],
[326859.17,326859.53,326860.25,326865.62,326874.56,326928.26,326999.89,327035.72,327071.56,327080.52,327089.49,327107.41,327125.34,327215.03,327304.79,327394.62,327754.68,328116.03,328660.75,329025.86,329576.72,329946.23,330317.65],
[378790.38,378790.82,378791.70,378798.32,378809.34,378875.51,378963.73,379007.84,379051.96,379062.99,379074.02,379096.08,379118.14,379228.44,379338.76,379449.09,379890.59,380332.43,380996.02,381439.08,382104.85,382549.59,382995.11],
[430721.40,430721.90,430722.89,430730.37,430742.84,430817.62,430917.31,430967.16,431017.00,431029.46,431041.92,431066.84,431091.76,431216.34,431340.91,431465.47,431963.59,432461.58,433208.44,433706.35,434453.36,434951.52,435449.87],
[482652.31,482652.84,482653.92,482662.02,482675.51,482756.43,482864.32,482918.25,482972.18,482985.66,482999.14,483026.11,483053.07,483187.86,483322.62,483457.36,483996.06,484534.39,485341.31,485878.91,486684.92,487222.04,487759.02],
[534583.13,534583.70,534584.84,534593.38,534607.62,534693.03,534806.89,534863.81,534920.73,534934.96,534949.19,534977.64,535006.10,535148.34,535290.56,535432.75,536001.18,536569.15,537420.29,537987.22,538836.94,539403.00,539968.77],
[586513.91,586514.50,586515.69,586524.56,586539.35,586628.09,586746.38,586805.52,586864.65,586879.44,586894.22,586923.78,586953.34,587101.13,587248.89,587396.62,587987.19,588577.26,589461.47,590050.39,590932.96,591520.84,592108.34],
[638444.65,638445.25,638446.47,638455.60,638470.80,638562.03,638683.65,638744.46,638805.26,638820.45,638835.65,638866.05,638896.44,639048.40,639200.31,639352.20,639959.43,640566.17,641475.38,642080.94,642988.46,643592.94,644197.01],
[690375.35,690375.98,690377.22,690386.53,690402.05,690495.16,690619.30,690681.36,690743.42,690758.93,690774.44,690805.47,690836.49,690991.59,691146.66,691301.69,691921.55,692540.93,693469.14,694087.39,695013.95,695631.14,696247.92],
[742306.04,742306.67,742307.93,742317.39,742333.14,742427.69,742553.72,742616.74,742679.74,742695.49,742711.24,742742.75,742774.25,742931.73,743089.18,743246.61,743876.04,744505.04,745447.73,746075.67,747016.81,747643.74,748270.28],
[794236.71,794237.35,794238.62,794248.19,794264.12,794359.74,794487.21,794550.95,794614.67,794630.60,794646.53,794678.40,794710.26,794869.54,795028.80,795188.03,795824.72,796461.01,797414.72,798050.04,799002.32,799636.71,800270.73],
[846167.37,846168.01,846169.30,846178.94,846195.01,846291.43,846419.98,846484.25,846548.52,846564.59,846580.65,846612.78,846644.91,846805.55,846966.16,847126.75,847768.90,848410.68,849372.71,850013.62,850974.35,851614.41,852254.14],
[898098.02,898098.66,898099.96,898109.66,898125.83,898222.85,898352.19,898416.86,898481.52,898497.69,898513.85,898546.18,898578.51,898740.15,898901.76,899063.36,899709.54,900355.41,901323.63,901968.73,902935.79,903580.11,904224.13],
[950028.66,950029.31,950030.61,950040.35,950056.59,950154.04,950283.95,950348.91,950413.86,950430.09,950446.33,950478.80,950511.28,950673.63,950835.97,950998.29,951647.40,952296.23,953268.96,953917.11,954888.82,955536.28,956183.47],
[1001959.29,1001959.94,1001961.25,1001971.02,1001987.31,1002085.05,1002215.36,1002280.51,1002345.66,1002361.94,1002378.23,1002410.80,1002443.37,1002606.23,1002769.06,1002931.88,1003583.02,1004233.92,1005209.82,1005860.12,1006835.12,1007484.81,1008134.27],
[1053889.92,1053890.57,1053891.88,1053901.67,1053917.99,1054015.92,1054146.48,1054211.75,1054277.03,1054293.35,1054309.66,1054342.30,1054374.93,1054538.10,1054701.26,1054864.41,1055516.85,1056169.10,1057147.07,1057798.80,1058775.99,1059427.20,1060078.19],
[1105820.54,1105821.20,1105822.50,1105832.31,1105848.65,1105946.67,1106077.36,1106142.71,1106208.05,1106224.39,1106240.72,1106273.39,1106306.06,1106469.41,1106632.74,1106796.06,1107449.23,1108102.24,1109081.41,1109733.97,1110712.48,1111364.59,1112016.51],
[1157751.17,1157751.82,1157753.13,1157762.93,1157779.28,1157877.33,1158008.07,1158073.43,1158138.79,1158155.14,1158171.48,1158204.16,1158236.84,1158400.24,1158563.63,1158727.01,1159380.44,1160033.72,1161013.37,1161666.29,1162645.37,1163297.91,1163950.28],
[1209681.78,1209682.44,1209683.75,1209693.55,1209709.89,1209807.91,1209938.61,1210003.96,1210069.31,1210085.64,1210101.98,1210134.65,1210167.33,1210330.68,1210494.03,1210657.38,1211310.68,1211963.85,1212943.40,1213596.27,1214575.34,1215227.90,1215880.32],
[1261612.40,1261613.05,1261614.36,1261624.16,1261640.48,1261738.43,1261869.04,1261934.33,1261999.63,1262015.96,1262032.28,1262064.93,1262097.58,1262260.82,1262424.05,1262587.28,1263240.13,1263892.88,1264871.83,1265524.34,1266502.90,1267155.15,1267807.29],
[1313543.02,1313543.67,1313544.97,1313554.76,1313571.06,1313668.90,1313799.36,1313864.58,1313929.80,1313946.11,1313962.42,1313995.03,1314027.64,1314190.69,1314353.74,1314516.79,1315168.93,1315821.00,1316798.95,1317450.83,1318428.49,1319080.16,1319731.74],
[1365473.63,1365474.28,1365475.58,1365485.35,1365501.64,1365599.33,1365729.59,1365794.72,1365859.85,1365876.13,1365892.41,1365924.98,1365957.54,1366120.36,1366283.17,1366445.99,1367097.20,1367748.36,1368725.00,1369376.03,1370352.44,1371003.31,1371654.10],
[1417404.24,1417404.89,1417406.19,1417415.94,1417432.20,1417529.72,1417659.76,1417724.77,1417789.79,1417806.04,1417822.30,1417854.81,1417887.31,1418049.85,1418212.39,1418374.92,1419025.04,1419675.12,1420650.17,1421300.16,1422275.05,1422924.93,1423574.75],
[1469334.85,1469335.50,1469336.80,1469346.53,1469362.76,1469460.09,1469589.87,1469654.76,1469719.65,1469735.87,1469752.09,1469784.54,1469816.98,1469979.20,1470141.42,1470303.64,1470952.51,1471601.37,1472574.61,1473223.41,1474196.56,1474845.29,1475493.99],
[1521265.46,1521266.11,1521267.41,1521277.12,1521293.31,1521390.43,1521519.93,1521584.68,1521649.43,1521665.62,1521681.81,1521714.18,1521746.56,1521908.44,1522070.31,1522232.19,1522879.70,1523527.20,1524498.45,1525145.94,1526117.16,1526764.62,1527412.06],
[1573196.07,1573196.72,1573198.01,1573207.70,1573223.85,1573320.75,1573449.96,1573514.56,1573579.16,1573595.31,1573611.46,1573643.77,1573676.07,1573837.57,1573999.08,1574160.59,1574806.64,1575452.70,1576421.81,1577067.89,1578037.01,1578683.10,1579329.18],
[1625126.68,1625127.33,1625128.62,1625138.28,1625154.39,1625251.06,1625379.95,1625444.40,1625508.85,1625524.96,1625541.07,1625573.29,1625605.52,1625766.64,1625927.76,1626088.88,1626733.39,1627377.92,1628344.77,1628989.35,1629956.26,1630600.89,1631245.52],
[1677057.29,1677057.93,1677059.22,1677068.86,1677084.93,1677181.36,1677309.92,1677374.21,1677438.49,1677454.56,1677470.64,1677502.78,1677534.92,1677695.64,1677856.35,1678017.08,1678659.98,1679302.93,1680267.40,1680910.43,1681875.02,1682518.11,1683161.23],
[1728987.90,1728988.54,1728989.82,1728999.44,1729015.47,1729111.64,1729239.88,1729303.99,1729368.11,1729384.14,1729400.17,1729432.23,1729464.29,1729624.59,1729784.89,1729945.20,1730586.45,1731227.75,1732189.79,1732831.20,1733793.38,1734434.89,1735076.43],
[1780918.51,1780919.15,1780920.43,1780930.02,1780946.00,1781041.92,1781169.81,1781233.76,1781297.70,1781313.69,1781329.68,1781361.65,1781393.63,1781553.50,1781713.38,1781873.26,1782512.82,1783152.44,1784111.97,1784751.72,1785711.44,1786351.31,1786991.22],
[1832849.11,1832849.75,1832851.03,1832860.59,1832876.54,1832972.19,1833099.74,1833163.51,1833227.28,1833243.22,1833259.17,1833291.06,1833322.94,1833482.38,1833641.83,1833801.27,1834439.11,1835077.01,1836033.99,1836672.05,1837629.25,1838267.45,1838905.71]])

#heat capacity ratio (=Cp,Cv)
gamma_PT = np.array(\
[[1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01,1.01,1.01,1.02,1.02,1.03],
[1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01,1.01,1.01,1.02,1.03,1.04,1.05],
[1.70,1.71,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01,1.01,1.01,1.02,1.03,1.04,1.05,1.07],
[1.69,1.70,1.70,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01,1.01,1.01,1.01,1.01,1.01,1.02,1.03,1.05,1.07,1.16,1.01,1.00],
[1.68,1.68,1.68,1.70,1.75,1.26,1.24,1.24,1.23,1.23,1.23,1.23,1.22,1.21,1.20,1.20,1.17,1.15,1.13,1.12,1.11,1.10,1.09],
[1.67,1.67,1.68,1.69,1.71,1.89,2.46,1.76,1.72,1.71,1.70,1.69,1.67,1.61,1.56,1.52,1.41,1.35,1.29,1.26,1.23,1.21,1.20],
[1.67,1.67,1.67,1.68,1.70,1.80,1.99,2.13,2.33,2.40,2.47,2.64,2.87,6.15,3.27,2.64,1.92,1.70,1.54,1.48,1.41,1.38,1.35],
[1.67,1.67,1.67,1.68,1.69,1.75,1.86,1.93,2.00,2.02,2.05,2.09,2.15,2.50,3.16,4.58,3.73,2.40,1.92,1.77,1.64,1.58,1.54],
[1.67,1.67,1.67,1.67,1.68,1.73,1.80,1.84,1.88,1.90,1.91,1.93,1.96,2.11,2.30,2.56,3.95,3.56,2.48,2.17,1.92,1.81,1.74],
[1.67,1.67,1.67,1.67,1.68,1.71,1.76,1.79,1.82,1.83,1.84,1.85,1.87,1.96,2.06,2.18,2.81,3.23,2.94,2.58,2.22,2.06,1.95],
[1.67,1.67,1.67,1.67,1.68,1.70,1.74,1.76,1.78,1.79,1.79,1.81,1.82,1.88,1.94,2.01,2.36,2.70,2.84,2.73,2.45,2.29,2.15],
[1.67,1.67,1.67,1.67,1.67,1.70,1.73,1.74,1.76,1.76,1.77,1.77,1.78,1.83,1.87,1.92,2.15,2.38,2.59,2.61,2.51,2.40,2.29],
[1.67,1.67,1.67,1.67,1.67,1.67,1.68,1.68,1.69,1.69,1.69,1.69,1.69,1.70,1.71,1.71,1.74,1.77,1.82,1.84,1.88,1.89,1.91],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.68,1.68,1.68,1.68,1.69,1.71,1.72,1.73,1.74,1.75,1.76],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.68,1.68,1.69,1.70,1.70,1.71,1.71],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.68,1.68,1.68,1.68,1.69,1.69],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.68,1.68,1.68],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66],
[1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.67,1.66,1.66]])

#Viscosity
mu_PT = np.array(\
[[1.35E-06,1.35E-06,1.35E-06,1.36E-06,1.36E-06,1.37E-06,1.38E-06,1.39E-06,1.40E-06,1.40E-06,1.40E-06,1.41E-06,1.41E-06,1.43E-06,1.45E-06,1.46E-06,1.53E-06,1.60E-06,1.69E-06,1.76E-06,1.85E-06,1.92E-06,2.00E-06],
[1.32E-06,1.32E-06,1.32E-06,1.32E-06,1.32E-06,1.33E-06,1.35E-06,1.36E-06,1.37E-06,1.37E-06,1.37E-06,1.38E-06,1.38E-06,1.40E-06,1.42E-06,1.44E-06,1.52E-06,1.60E-06,1.72E-06,1.80E-06,1.94E-06,2.06E-06,2.19E-06],
[4.61E-07,4.62E-07,1.30E-06,1.30E-06,1.30E-06,1.32E-06,1.34E-06,1.35E-06,1.36E-06,1.36E-06,1.36E-06,1.37E-06,1.37E-06,1.39E-06,1.42E-06,1.44E-06,1.54E-06,1.64E-06,1.81E-06,1.94E-06,2.18E-06,2.39E-06,2.66E-06],
[5.11E-07,5.12E-07,5.13E-07,1.49E-06,1.49E-06,1.51E-06,1.54E-06,1.55E-06,1.57E-06,1.57E-06,1.58E-06,1.58E-06,1.59E-06,1.63E-06,1.67E-06,1.72E-06,1.92E-06,2.17E-06,2.68E-06,3.16E-06,4.28E-06,4.90E-06,5.41E-06],
[7.64E-07,7.64E-07,7.65E-07,7.70E-07,7.79E-07,3.58E-06,3.62E-06,3.63E-06,3.65E-06,3.66E-06,3.66E-06,3.67E-06,3.68E-06,3.73E-06,3.79E-06,3.84E-06,4.08E-06,4.32E-06,4.69E-06,4.93E-06,5.28E-06,5.50E-06,5.71E-06],
[1.02E-06,1.02E-06,1.02E-06,1.03E-06,1.03E-06,1.07E-06,1.15E-06,3.33E-06,3.38E-06,3.39E-06,3.40E-06,3.42E-06,3.45E-06,3.56E-06,3.66E-06,3.76E-06,4.13E-06,4.47E-06,4.96E-06,5.28E-06,5.75E-06,6.06E-06,6.38E-06],
[1.25E-06,1.25E-06,1.25E-06,1.26E-06,1.26E-06,1.30E-06,1.36E-06,1.39E-06,1.43E-06,1.44E-06,1.45E-06,1.47E-06,1.49E-06,2.68E-06,2.95E-06,3.11E-06,3.57E-06,3.92E-06,4.38E-06,4.66E-06,5.07E-06,5.34E-06,5.60E-06],
[1.46E-06,1.46E-06,1.46E-06,1.47E-06,1.47E-06,1.51E-06,1.55E-06,1.58E-06,1.61E-06,1.61E-06,1.62E-06,1.64E-06,1.65E-06,1.73E-06,1.84E-06,1.99E-06,2.94E-06,3.42E-06,3.90E-06,4.17E-06,4.55E-06,4.79E-06,5.02E-06],
[1.66E-06,1.66E-06,1.66E-06,1.66E-06,1.67E-06,1.70E-06,1.74E-06,1.76E-06,1.78E-06,1.79E-06,1.80E-06,1.81E-06,1.82E-06,1.88E-06,1.95E-06,2.03E-06,2.44E-06,2.95E-06,3.51E-06,3.79E-06,4.15E-06,4.37E-06,4.59E-06],
[1.84E-06,1.84E-06,1.84E-06,1.85E-06,1.85E-06,1.88E-06,1.92E-06,1.93E-06,1.95E-06,1.96E-06,1.96E-06,1.97E-06,1.98E-06,2.04E-06,2.09E-06,2.15E-06,2.41E-06,2.73E-06,3.21E-06,3.50E-06,3.86E-06,4.07E-06,4.27E-06],
[2.02E-06,2.02E-06,2.02E-06,2.02E-06,2.02E-06,2.05E-06,2.08E-06,2.10E-06,2.12E-06,2.12E-06,2.13E-06,2.13E-06,2.14E-06,2.19E-06,2.23E-06,2.28E-06,2.49E-06,2.72E-06,3.08E-06,3.32E-06,3.66E-06,3.86E-06,4.05E-06],
[2.18E-06,2.18E-06,2.18E-06,2.19E-06,2.19E-06,2.21E-06,2.24E-06,2.26E-06,2.27E-06,2.28E-06,2.28E-06,2.29E-06,2.30E-06,2.34E-06,2.38E-06,2.42E-06,2.59E-06,2.78E-06,3.07E-06,3.26E-06,3.55E-06,3.73E-06,3.90E-06],
[3.54E-06,3.54E-06,3.54E-06,3.54E-06,3.55E-06,3.56E-06,3.57E-06,3.58E-06,3.59E-06,3.59E-06,3.59E-06,3.60E-06,3.60E-06,3.62E-06,3.64E-06,3.66E-06,3.74E-06,3.82E-06,3.93E-06,4.01E-06,4.12E-06,4.20E-06,4.27E-06],
[4.61E-06,4.61E-06,4.61E-06,4.61E-06,4.61E-06,4.62E-06,4.63E-06,4.63E-06,4.64E-06,4.64E-06,4.64E-06,4.64E-06,4.65E-06,4.66E-06,4.68E-06,4.69E-06,4.74E-06,4.80E-06,4.88E-06,4.93E-06,5.01E-06,5.06E-06,5.12E-06],
[5.52E-06,5.52E-06,5.52E-06,5.52E-06,5.52E-06,5.53E-06,5.54E-06,5.54E-06,5.55E-06,5.55E-06,5.55E-06,5.55E-06,5.55E-06,5.56E-06,5.58E-06,5.59E-06,5.63E-06,5.67E-06,5.74E-06,5.78E-06,5.85E-06,5.89E-06,5.93E-06],
[6.34E-06,6.34E-06,6.34E-06,6.34E-06,6.34E-06,6.35E-06,6.36E-06,6.36E-06,6.36E-06,6.37E-06,6.37E-06,6.37E-06,6.37E-06,6.38E-06,6.39E-06,6.40E-06,6.44E-06,6.47E-06,6.53E-06,6.57E-06,6.62E-06,6.66E-06,6.70E-06],
[7.10E-06,7.10E-06,7.10E-06,7.10E-06,7.10E-06,7.11E-06,7.11E-06,7.12E-06,7.12E-06,7.12E-06,7.12E-06,7.12E-06,7.13E-06,7.13E-06,7.14E-06,7.15E-06,7.18E-06,7.22E-06,7.27E-06,7.30E-06,7.35E-06,7.39E-06,7.42E-06],
[7.81E-06,7.81E-06,7.81E-06,7.81E-06,7.81E-06,7.82E-06,7.82E-06,7.83E-06,7.83E-06,7.83E-06,7.83E-06,7.83E-06,7.84E-06,7.84E-06,7.85E-06,7.86E-06,7.89E-06,7.92E-06,7.97E-06,8.00E-06,8.04E-06,8.07E-06,8.10E-06],
[8.49E-06,8.49E-06,8.49E-06,8.49E-06,8.49E-06,8.49E-06,8.50E-06,8.50E-06,8.51E-06,8.51E-06,8.51E-06,8.51E-06,8.51E-06,8.52E-06,8.53E-06,8.53E-06,8.56E-06,8.59E-06,8.63E-06,8.66E-06,8.70E-06,8.73E-06,8.76E-06],
[9.14E-06,9.14E-06,9.14E-06,9.14E-06,9.14E-06,9.14E-06,9.15E-06,9.15E-06,9.15E-06,9.16E-06,9.16E-06,9.16E-06,9.16E-06,9.17E-06,9.17E-06,9.18E-06,9.21E-06,9.23E-06,9.27E-06,9.30E-06,9.34E-06,9.37E-06,9.39E-06],
[9.77E-06,9.77E-06,9.77E-06,9.77E-06,9.77E-06,9.77E-06,9.78E-06,9.78E-06,9.78E-06,9.78E-06,9.78E-06,9.78E-06,9.78E-06,9.79E-06,9.80E-06,9.80E-06,9.83E-06,9.85E-06,9.89E-06,9.92E-06,9.95E-06,9.98E-06,1.00E-05],
[1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.02E-05,1.03E-05,1.03E-05,1.03E-05,1.04E-05,1.04E-05,1.04E-05],
[1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.08E-05,1.09E-05,1.09E-05,1.09E-05,1.10E-05,1.10E-05,1.10E-05],
[1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.14E-05,1.15E-05,1.15E-05,1.15E-05,1.15E-05,1.16E-05],
[1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.19E-05,1.20E-05,1.20E-05,1.20E-05,1.20E-05,1.20E-05,1.20E-05,1.21E-05,1.21E-05,1.21E-05,1.21E-05],
[1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.25E-05,1.26E-05,1.26E-05,1.26E-05,1.26E-05,1.27E-05,1.27E-05],
[1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.30E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.31E-05,1.32E-05,1.32E-05,1.32E-05],
[1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.36E-05,1.37E-05,1.37E-05,1.37E-05,1.37E-05,1.37E-05],
[1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.41E-05,1.42E-05,1.42E-05,1.42E-05,1.42E-05,1.42E-05,1.43E-05],
[1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.46E-05,1.47E-05,1.47E-05,1.47E-05,1.47E-05,1.47E-05,1.48E-05,1.48E-05],
[1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.51E-05,1.52E-05,1.52E-05,1.52E-05,1.52E-05,1.52E-05,1.52E-05,1.52E-05,1.53E-05,1.53E-05],
[1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.56E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.57E-05,1.58E-05,1.58E-05],
[1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.61E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.62E-05,1.63E-05],
[1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.66E-05,1.67E-05,1.67E-05,1.67E-05,1.67E-05,1.67E-05,1.67E-05,1.67E-05],
[1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.71E-05,1.72E-05,1.72E-05,1.72E-05,1.72E-05,1.72E-05,1.72E-05],
[1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.76E-05,1.77E-05,1.77E-05,1.77E-05,1.77E-05,1.77E-05],
[1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.81E-05,1.82E-05,1.82E-05],
[1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.85E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05,1.86E-05],
[1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.90E-05,1.91E-05,1.91E-05,1.91E-05,1.91E-05,1.91E-05],
[1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05,1.95E-05],
[1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,1.99E-05,2.00E-05,2.00E-05,2.00E-05,2.00E-05,2.00E-05,2.00E-05],
[2.04E-05,2.04E-05,2.04E-05,2.04E-05,2.04E-05,2.03E-05,2.03E-05,2.02E-05,2.02E-05,2.02E-05,2.02E-05,2.02E-05,2.02E-05,2.01E-05,2.00E-05,1.99E-05,1.97E-05,1.94E-05,1.90E-05,1.87E-05,1.83E-05,1.80E-05,1.77E-05],
[2.08E-05,2.08E-05,2.08E-05,2.08E-05,2.08E-05,2.08E-05,2.07E-05,2.07E-05,2.07E-05,2.07E-05,2.06E-05,2.06E-05,2.06E-05,2.05E-05,2.05E-05,2.04E-05,2.01E-05,1.99E-05,1.95E-05,1.92E-05,1.88E-05,1.85E-05,1.83E-05],
[2.13E-05,2.13E-05,2.13E-05,2.13E-05,2.13E-05,2.12E-05,2.12E-05,2.11E-05,2.11E-05,2.11E-05,2.11E-05,2.11E-05,2.11E-05,2.10E-05,2.09E-05,2.09E-05,2.06E-05,2.03E-05,1.99E-05,1.97E-05,1.93E-05,1.90E-05,1.88E-05],
[2.17E-05,2.17E-05,2.17E-05,2.17E-05,2.17E-05,2.17E-05,2.16E-05,2.16E-05,2.16E-05,2.15E-05,2.15E-05,2.15E-05,2.15E-05,2.14E-05,2.14E-05,2.13E-05,2.11E-05,2.08E-05,2.04E-05,2.02E-05,1.98E-05,1.95E-05,1.93E-05],
[2.21E-05,2.21E-05,2.21E-05,2.21E-05,2.21E-05,2.21E-05,2.20E-05,2.20E-05,2.20E-05,2.20E-05,2.20E-05,2.20E-05,2.20E-05,2.19E-05,2.18E-05,2.18E-05,2.15E-05,2.13E-05,2.09E-05,2.06E-05,2.03E-05,2.00E-05,1.98E-05]])


interp_P_T_hPT = interp2d(P,T,h_PT)
interp_P_T_DPT = interp2d(P,T,D_PT)
interp_P_T_mu = interp2d(P,T,mu_PT)

