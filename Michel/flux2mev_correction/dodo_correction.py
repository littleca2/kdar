"""
//Usage
//double ReferenceFlux8MeV() return: double (reference flux conversion factor for 8 MeV)
//int DatasetNumberCheck(int RUN, int SUBRUN) return: int (0-99, error: -1)
//double TimeCorrrection(int TimePeriod) return: double (correction factor, error: -1)
//double SpatialCorrrection(int R2, int Z) return: double (correction factor, error: -1)
"""
def ReferenceFlux8MeV():
    flux = 7420.31 #PSD volume, Flux for 8MeV
    #flux = 7766.54 //for 2022 (this difference is included in TimeCorrrection part)
    #flux = 7366.86 //for 2023 (this difference is included in TimeCorrection part)

    return flux


def DatasetNumberCheck(RUN, SUBRUN):
    #delayed Run/Subrun ID will be used
    num = -1
    tmpRUN = RUN*100000+SUBRUN

    if(tmpRUN < 160002532):
         num = 0
    elif(tmpRUN >= 160002532 and tmpRUN < 160902818):
         num = 1
    elif(tmpRUN >= 160902818 and tmpRUN < 162501792):
         num = 2
    elif(tmpRUN >= 162501792 and tmpRUN < 162506792):
         num = 3
    elif(tmpRUN >= 162506792 and tmpRUN < 162600916):
         num = 4
    elif(tmpRUN >= 162600916 and tmpRUN < 163400811):
         num = 5
    elif(tmpRUN >= 163400811 and tmpRUN < 167200638):
         num = 6
    elif(tmpRUN >= 167200638 and tmpRUN < 167800545):
         num = 7
    elif(tmpRUN >= 167800545 and tmpRUN < 167805545):
         num = 8
    elif(tmpRUN >= 167805545 and tmpRUN < 168401608):
         num = 9
    elif(tmpRUN >= 168401608 and tmpRUN < 168406608):
         num = 10
    elif(tmpRUN >= 168406608 and tmpRUN < 168900010):
         num = 11
    elif(tmpRUN >= 168900010 and tmpRUN < 168905010):
         num = 12
    elif(tmpRUN >= 168905010 and tmpRUN < 168910010):
         num = 13
    elif(tmpRUN >= 168910010 and tmpRUN < 169201720):
         num = 14
    elif(tmpRUN >= 169201720 and tmpRUN < 171701997):
         num = 15
    elif(tmpRUN >= 171701997 and tmpRUN < 171801105):
         num = 16
    elif(tmpRUN >= 171801105 and tmpRUN < 171806105):
         num = 17
    elif(tmpRUN >= 171806105 and tmpRUN < 172201282):
         num = 18
    elif(tmpRUN >= 172201282 and tmpRUN < 172500331):
         num = 19
    elif(tmpRUN >= 172500331 and tmpRUN < 172601638):
         num = 20
    elif(tmpRUN >= 172601638 and tmpRUN < 172903435):
         num = 21
    elif(tmpRUN >= 172903435 and tmpRUN < 174903245):
         num = 22
    elif(tmpRUN >= 174903245 and tmpRUN < 176900682):
         num = 23
    elif(tmpRUN >= 176900682 and tmpRUN < 178501725):
         num = 24
    elif(tmpRUN >= 178501725 and tmpRUN < 178704458):
         num = 25
    elif(tmpRUN >= 178704458 and tmpRUN < 178804148):
         num = 26
    elif(tmpRUN >= 178804148 and tmpRUN < 179003501):
         num = 27
    elif(tmpRUN >= 179003501 and tmpRUN < 181400826):
         num = 28
    elif(tmpRUN >= 181400826 and tmpRUN < 181504483):
         num = 29
    elif(tmpRUN >= 181504483 and tmpRUN < 181801197):
         num = 30
    elif(tmpRUN >= 181801197 and tmpRUN < 182601544):
         num = 31
    elif(tmpRUN >= 182601544 and tmpRUN < 182802561):
         num = 32
    elif(tmpRUN >= 182802561 and tmpRUN < 183002420):
         num = 33
    elif(tmpRUN >= 183002420 and tmpRUN < 183101359):
         num = 34
    elif(tmpRUN >= 183101359 and tmpRUN < 183700994):
         num = 35
    elif(tmpRUN >= 183700994 and tmpRUN < 183902723):
         num = 36
    elif(tmpRUN >= 183902723 and tmpRUN < 183907723):
         num = 37
    elif(tmpRUN >= 183907723 and tmpRUN < 184304381):
         num = 38
    elif(tmpRUN >= 184304381 and tmpRUN < 185003641):
         num = 39
    elif(tmpRUN >= 185003641 and tmpRUN < 185008641):
         num = 40
    elif(tmpRUN >= 185008641 and tmpRUN < 185013641):
         num = 41
    elif(tmpRUN >= 185013641 and tmpRUN < 186400358):
         num = 42
    elif(tmpRUN >= 186400358 and tmpRUN < 186700240):
         num = 43
    elif(tmpRUN >= 186700240 and tmpRUN < 187400303):
         num = 44
    elif(tmpRUN >= 187400303 and tmpRUN < 187802861):
         num = 45
    elif(tmpRUN >= 187802861 and tmpRUN < 187807861):
         num = 46
    elif(tmpRUN >= 187807861 and tmpRUN < 187901356):
         num = 47
    elif(tmpRUN >= 187901356 and tmpRUN < 190902370):
         num = 48
    elif(tmpRUN >= 190902370 and tmpRUN < 190907370):
         num = 49
    elif(tmpRUN >= 190907370 and tmpRUN < 191004199):
         num = 50
    elif(tmpRUN >= 191004199 and tmpRUN < 191401909):
         num = 51
    elif(tmpRUN >= 191401909 and tmpRUN < 191603447):
         num = 52
    elif(tmpRUN >= 191603447 and tmpRUN < 191608447):
         num = 53
    elif(tmpRUN >= 191608447 and tmpRUN < 191702034):
         num = 54
    elif(tmpRUN >= 191702034 and tmpRUN < 192200152):
         num = 55
    elif(tmpRUN >= 192200152 and tmpRUN < 192303167):
         num = 56
    elif(tmpRUN >= 192303167 and tmpRUN < 193000055):
         num = 57
    elif(tmpRUN >= 193000055 and tmpRUN < 193601809):
         num = 58
    elif(tmpRUN >= 193601809 and tmpRUN < 194003595):
         num = 59
    elif(tmpRUN >= 194003595 and tmpRUN < 194401518):
         num = 60
    elif(tmpRUN >= 194401518 and tmpRUN < 194406518):
         num = 61
    elif(tmpRUN >= 194406518 and tmpRUN < 195000194):
         num = 62
    elif(tmpRUN >= 195000194 and tmpRUN < 195201349):
         num = 63
    elif(tmpRUN >= 195201349 and tmpRUN < 198099999):
         num = 64 #2021: -run 1980
    elif(tmpRUN >= 198100000 and tmpRUN < 209001219):
         num = 65 #2022: run 1981- (2070?-)
    elif(tmpRUN >= 209001219 and tmpRUN < 209006219):
         num = 66
    elif(tmpRUN >= 209006219 and tmpRUN < 209503996):
         num = 67
    elif(tmpRUN >= 209503996 and tmpRUN < 209604315):
         num = 68
    elif(tmpRUN >= 209604315 and tmpRUN < 209802190):
         num = 69
    elif(tmpRUN >= 209802190 and tmpRUN < 210004954):
         num = 70
    elif(tmpRUN >= 210004954 and tmpRUN < 211302031):
         num = 71
    elif(tmpRUN >= 211302031 and tmpRUN < 211307031):
         num = 72
    elif(tmpRUN >= 211307031 and tmpRUN < 212200485):
         num = 73
    elif(tmpRUN >= 212200485 and tmpRUN < 212205485):
         num = 74
    elif(tmpRUN >= 212205485 and tmpRUN < 212602873):
         num = 75
    elif(tmpRUN >= 212602873 and tmpRUN < 213401750):
         num = 76
    elif(tmpRUN >= 213401750 and tmpRUN < 213801039):
         num = 77
    elif(tmpRUN >= 213801039 and tmpRUN < 213806039):
         num = 78
    elif(tmpRUN >= 213806039 and tmpRUN < 214502387):
         num = 79
    elif(tmpRUN >= 214502387 and tmpRUN < 214602600):
         num = 80
    elif(tmpRUN >= 214602600 and tmpRUN < 214703772):
         num = 81
    elif(tmpRUN >= 214703772 and tmpRUN < 215802016):
         num = 82
    elif(tmpRUN >= 215802016 and tmpRUN < 216803160):
         num = 83
    elif(tmpRUN >= 216803160 and tmpRUN < 216900406):
         num = 84
    elif(tmpRUN >= 216900406 and tmpRUN < 217000858):
         num = 85
    elif(tmpRUN >= 217000858 and tmpRUN < 217403444):
         num = 86
    elif(tmpRUN >= 217403444 and tmpRUN < 217408444):
         num = 87
    elif(tmpRUN >= 217408444 and tmpRUN < 217413444):
         num = 88
    elif(tmpRUN >= 217413444 and tmpRUN < 218000393):
         num = 89
    elif(tmpRUN >= 218000393 and tmpRUN < 218104400):
         num = 90
    elif(tmpRUN >= 218104400 and tmpRUN < 218109400):
         num = 91
    elif(tmpRUN >= 218109400 and tmpRUN < 218114400):
         num = 92
    elif(tmpRUN >= 218114400 and tmpRUN < 218903593):
         num = 93
    elif(tmpRUN >= 218903593 and tmpRUN < 218908593):
         num = 94
    elif(tmpRUN >= 218908593 and tmpRUN < 218913593):
         num = 95
    elif(tmpRUN >= 218913593 and tmpRUN < 219501487):
         num = 96
    elif(tmpRUN >= 219501487 and tmpRUN < 219604719):
         num = 97
    elif(tmpRUN >= 219604719 and tmpRUN < 219900370):
         num = 98
    elif(tmpRUN >= 219900370 and tmpRUN < 221099999):
         num = 99 #2022: -run 2210
    #2023 sterile:run 2248- (kicker/self: run 2226?-)
	# This is not up to date
#    elif(tmpRUN >= 224800000 and tmpRUN < 225902198):
#         num = 100
#    elif(tmpRUN >= 225902198 and tmpRUN < 225907198):
#         num = 101
#    elif(tmpRUN >= 225907198 and tmpRUN < 226000538):
#         num = 102
#    elif(tmpRUN >= 226000538 and tmpRUN < 227200393):
#         num = 103
#    elif(tmpRUN >= 227200393 and tmpRUN < 227804386):
#         num = 104
#    elif(tmpRUN >= 227804386 and tmpRUN < 228300747):
#         num = 105
#    elif(tmpRUN >= 228300747 and tmpRUN < 228305747):
#        num = 106
#    elif(tmpRUN >= 228305747 and tmpRUN < 228310747):
#        num = 107
#    elif(tmpRUN >= 228310747 and tmpRUN < 229000392):
#        num = 108
#    elif(tmpRUN >= 229000392 and tmpRUN < 229005392):
#        num = 109
#    elif(tmpRUN >= 229005392 and tmpRUN < 229010392):
#        num = 110
#    elif(tmpRUN >= 229010392 and tmpRUN < 229399999):
#        num = 111
    
    return num

def TimeCorrection(TimePeriod):
    val = -1.0
    if(TimePeriod == 0):
        val = 1.00000
    elif(TimePeriod == 1):
        val = 1.00182
    elif(TimePeriod == 2):
        val = 1.00580
    elif(TimePeriod == 3):
        val = 1.00486
    elif(TimePeriod == 4):
        val = 1.00604
    elif(TimePeriod == 5):
        val = 1.00505
    elif(TimePeriod == 6):
        val = 1.00625
    elif(TimePeriod == 7):
        val = 1.00466
    elif(TimePeriod == 8):
        val = 1.00406
    elif(TimePeriod == 9):
        val = 1.00375
    elif(TimePeriod == 10):
        val = 1.00273
    elif(TimePeriod == 11):
        val = 1.00395
    elif(TimePeriod == 12):
        val = 1.00644
    elif(TimePeriod == 13):
        val = 1.00338
    elif(TimePeriod == 14):
        val = 1.00691
    elif(TimePeriod == 15):
        val = 1.00915
    elif(TimePeriod == 16):
        val = 1.01197
    elif(TimePeriod == 17):
        val = 1.01135
    elif(TimePeriod == 18):
        val = 1.01209
    elif(TimePeriod == 19):
        val = 1.01121
    elif(TimePeriod == 20):
        val = 1.01254
    elif(TimePeriod == 21):
        val = 1.01154
    elif(TimePeriod == 22):
        val = 1.01460
    elif(TimePeriod == 23):
        val = 1.01080
    elif(TimePeriod == 24):
        val = 1.00911
    elif(TimePeriod == 25):
        val = 1.01003
    elif(TimePeriod == 26):
        val = 1.01104
    elif(TimePeriod == 27):
        val = 1.01153
    elif(TimePeriod == 28):
        val = 1.01248
    elif(TimePeriod == 29):
        val = 1.00921
    elif(TimePeriod == 30):
        val = 1.00951
    elif(TimePeriod == 31):
        val = 1.01118
    elif(TimePeriod == 32):
        val = 1.01114
    elif(TimePeriod == 33):
        val = 1.01210
    elif(TimePeriod == 34):
        val = 1.01458
    elif(TimePeriod == 35):
        val = 1.01308
    elif(TimePeriod == 36):
        val = 1.01577
    elif(TimePeriod == 37):
        val = 1.01571
    elif(TimePeriod == 38):
        val = 1.02102
    elif(TimePeriod == 39):
        val = 1.02025
    elif(TimePeriod == 40):
        val = 1.02026
    elif(TimePeriod == 41):
        val = 1.02502
    elif(TimePeriod == 42):
        val = 1.02644
    elif(TimePeriod == 43):
        val = 1.02869
    elif(TimePeriod == 44):
        val = 1.02494
    elif(TimePeriod == 45):
        val = 1.02688
    elif(TimePeriod == 46):
        val = 1.02610
    elif(TimePeriod == 47):
        val = 1.02415
    elif(TimePeriod == 48):
        val = 1.02553
    elif(TimePeriod == 49):
        val = 1.02369
    elif(TimePeriod == 50):
        val = 1.02243
    elif(TimePeriod == 51):
        val = 1.02243
    elif(TimePeriod == 52):
        val = 1.02531
    elif(TimePeriod == 53):
        val = 1.02344
    elif(TimePeriod == 54):
        val = 1.02223
    elif(TimePeriod == 55):
        val = 1.02324
    elif(TimePeriod == 56):
        val = 1.02255
    elif(TimePeriod == 57):
        val = 1.02148
    elif(TimePeriod == 58):
        val = 1.02211
    elif(TimePeriod == 59):
        val = 1.01947
    elif(TimePeriod == 60):
        val = 1.01794
    elif(TimePeriod == 61):
        val = 1.01917
    elif(TimePeriod == 62):
        val = 1.01657
    elif(TimePeriod == 63):
        val = 1.01760
    elif(TimePeriod == 64):
        val = 1.01794		#2021
    elif(TimePeriod == 65):
        val = 1.00000		#2022
    elif(TimePeriod == 66):
        val = 0.99801
    elif(TimePeriod == 67):
        val = 0.99897
    elif(TimePeriod == 68):
        val = 0.99833
    elif(TimePeriod == 69):
        val = 0.99705
    elif(TimePeriod == 70):
        val = 0.99659
    elif(TimePeriod == 71):
        val = 0.99500
    elif(TimePeriod == 72):
        val = 0.99382
    elif(TimePeriod == 73):
        val = 0.99391
    elif(TimePeriod == 74):
        val = 0.99103
    elif(TimePeriod == 75):
        val = 0.99199
    elif(TimePeriod == 76):
        val = 0.99278
    elif(TimePeriod == 77):
        val = 0.99585
    elif(TimePeriod == 78):
        val = 0.99529
    elif(TimePeriod == 79):
        val = 0.99321
    elif(TimePeriod == 80):
        val = 0.99407
    elif(TimePeriod == 81):
        val = 0.99298
    elif(TimePeriod == 82):
        val = 0.99504
    elif(TimePeriod == 83):
        val = 0.99355
    elif(TimePeriod == 84):
        val = 0.99437
    elif(TimePeriod == 85):
        val = 0.99472
    elif(TimePeriod == 86):
        val = 0.99614
    elif(TimePeriod == 87):
        val = 0.99655
    elif(TimePeriod == 88):
        val = 0.99417
    elif(TimePeriod == 89):
        val = 0.99347
    elif(TimePeriod == 90):
        val = 0.99324
    elif(TimePeriod == 91):
        val = 0.99420
    elif(TimePeriod == 92):
        val = 0.99215
    elif(TimePeriod == 93):
        val = 0.99325
    elif(TimePeriod == 94):
        val = 0.99021
    elif(TimePeriod == 95):
        val = 0.98920
    elif(TimePeriod == 96):
        val = 0.98646
    elif(TimePeriod == 97):
        val = 0.98643
    elif(TimePeriod == 98):
        val = 0.98369
    elif(TimePeriod == 99):
        val = 0.98534		#2022
    elif(TimePeriod == 100):
        val = 1.00000		#2023
    elif(TimePeriod == 101):
        val = 0.99715
    elif(TimePeriod == 102):
        val = 0.99779
    elif(TimePeriod == 103):
        val = 1.00166
    elif(TimePeriod == 104):
        val = 1.00250
    elif(TimePeriod == 105):
        val = 1.00012
    elif(TimePeriod == 106):
        val = 0.99804
    elif(TimePeriod == 107):
        val = 0.99614
    elif(TimePeriod == 108):
        val = 0.99768
    elif(TimePeriod == 109):
        val = 0.99579
    elif(TimePeriod == 110):
        val = 0.99553
    elif(TimePeriod == 111):
        val = 0.99347//2023
        
    if(TimePeriod >= 65 and TimePeriod < 100):
        val *= (7766.54/7420.31)
    if(TimePeriod >= 100 and TimePeriod < 112):
        val *= (7366.86/7420.31)

    return val

"""
int R2Location(double R2m2)
{
    int num = -1
    if(R2m2 < 0.218) num = 0
    elif(R2m2 >= 0.218 and R2m2 < 0.436) num = 1
    elif(R2m2 >= 0.436 and R2m2 < 0.654) num = 2
    elif(R2m2 >= 0.654 and R2m2 < 0.872) num = 3
    elif(R2m2 >= 0.872 and R2m2 < 1.090) num = 4
    elif(R2m2 >= 1.090 and R2m2 < 1.308) num = 5
    elif(R2m2 >= 1.308 and R2m2 < 1.526) num = 6
    elif(R2m2 >= 1.526 and R2m2 < 1.744) num = 7
    elif(R2m2 >= 1.744 and R2m2 < 1.96) num = 8
    //elif(R2m2 >= 1.744 and R2m2 < 1.962) num = 8
    return num
}


int ZLocation(double Zm)
{
    int num = -1
    if(Zm >= -1.0 and Zm < -0.75) num = 0
    elif(Zm >= -0.75 and Zm < -0.50) num = 1
    elif(Zm >= -0.50 and Zm < -0.25) num = 2
    elif(Zm >= -0.25 and Zm < 0.0) num = 3
    elif(Zm >= 0.0 and Zm < 0.25) num = 4
    elif(Zm >= 0.25 and Zm < 0.50) num = 5
    elif(Zm >= 0.50 and Zm < 0.75) num = 6
    elif(Zm >= 0.75 and Zm < 1.00) num = 7
    return num
}


double SpatialCorrrection(int R2, int Z)
{
    double val = -1.0
    if(R2 == 0 and Z == 0) val = 1.00826
    elif(R2 == 0 and Z == 1) val = 1.01131
    elif(R2 == 0 and Z == 2) val = 1.01326
    elif(R2 == 0 and Z == 3) val = 1.01435
    elif(R2 == 0 and Z == 4) val = 1.01266
    elif(R2 == 0 and Z == 5) val = 1.00963
    elif(R2 == 0 and Z == 6) val = 1.00295
    elif(R2 == 0 and Z == 7) val = 0.988422
    elif(R2 == 1 and Z == 0) val = 1.00082
    elif(R2 == 1 and Z == 1) val = 1.0104
    elif(R2 == 1 and Z == 2) val = 1.01332
    elif(R2 == 1 and Z == 3) val = 1.01207
    elif(R2 == 1 and Z == 4) val = 1.01006
    elif(R2 == 1 and Z == 5) val = 1.00594
    elif(R2 == 1 and Z == 6) val = 1.00291
    elif(R2 == 1 and Z == 7) val = 0.988636
    elif(R2 == 2 and Z == 0) val = 1.0035
    elif(R2 == 2 and Z == 1) val = 1.00646
    elif(R2 == 2 and Z == 2) val = 1.00816
    elif(R2 == 2 and Z == 3) val = 1.00943
    elif(R2 == 2 and Z == 4) val = 1.00865
    elif(R2 == 2 and Z == 5) val = 1.00674
    elif(R2 == 2 and Z == 6) val = 1.00121
    elif(R2 == 2 and Z == 7) val = 0.986419
    elif(R2 == 3 and Z == 0) val = 1.00133
    elif(R2 == 3 and Z == 1) val = 1.00654
    elif(R2 == 3 and Z == 2) val = 1.00601
    elif(R2 == 3 and Z == 3) val = 1.00455
    elif(R2 == 3 and Z == 4) val = 1.00384
    elif(R2 == 3 and Z == 5) val = 1.00407
    elif(R2 == 3 and Z == 6) val = 1.00039
    elif(R2 == 3 and Z == 7) val = 0.990121
    elif(R2 == 4 and Z == 0) val = 1.00042
    elif(R2 == 4 and Z == 1) val = 1.00354
    elif(R2 == 4 and Z == 2) val = 1.00401
    elif(R2 == 4 and Z == 3) val = 1.00266
    elif(R2 == 4 and Z == 4) val = 1.00035
    elif(R2 == 4 and Z == 5) val = 0.9993
    elif(R2 == 4 and Z == 6) val = 0.99988
    elif(R2 == 4 and Z == 7) val = 0.990895
    elif(R2 == 5 and Z == 0) val = 0.994486
    elif(R2 == 5 and Z == 1) val = 1.00247
    elif(R2 == 5 and Z == 2) val = 0.999726
    elif(R2 == 5 and Z == 3) val = 0.996691
    elif(R2 == 5 and Z == 4) val = 0.9951
    elif(R2 == 5 and Z == 5) val = 0.995753
    elif(R2 == 5 and Z == 6) val = 0.996516
    elif(R2 == 5 and Z == 7) val = 0.98747
    elif(R2 == 6 and Z == 0) val = 0.995705
    elif(R2 == 6 and Z == 1) val = 1.00121
    elif(R2 == 6 and Z == 2) val = 0.990488
    elif(R2 == 6 and Z == 3) val = 0.993217
    elif(R2 == 6 and Z == 4) val = 0.987989
    elif(R2 == 6 and Z == 5) val = 0.992345
    elif(R2 == 6 and Z == 6) val = 0.989659
    elif(R2 == 6 and Z == 7) val = 0.980401
    elif(R2 == 7 and Z == 0) val = 1.00092
    elif(R2 == 7 and Z == 1) val = 0.996873
    elif(R2 == 7 and Z == 2) val = 0.992227
    elif(R2 == 7 and Z == 3) val = 0.985249
    elif(R2 == 7 and Z == 4) val = 0.991992
    elif(R2 == 7 and Z == 5) val = 0.99511
    elif(R2 == 7 and Z == 6) val = 0.987255
    elif(R2 == 7 and Z == 7) val = 0.982257
    elif(R2 == 8 and Z == 0) val = 1.00062
    elif(R2 == 8 and Z == 1) val = 0.997627
    elif(R2 == 8 and Z == 2) val = 0.991336
    elif(R2 == 8 and Z == 3) val = 0.977706
    elif(R2 == 8 and Z == 4) val = 0.988281
    elif(R2 == 8 and Z == 5) val = 0.975438
    elif(R2 == 8 and Z == 6) val = 0.988966
    elif(R2 == 8 and Z == 7) val = 0.978012

    return val
}

"""

