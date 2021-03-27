#algoritm genetic pentru determinarea maximului unei funcţii pozitive pe un domeniu dat
import random
from math import log
from bitarray import bitarray, util
from random import randint, choice, uniform
import matplotlib.pyplot as plt

f = open("input.txt", "r")
g = open("output.txt", "w")

#citesc datele de intrare din fisier
nr_populatie = int(f.readline())
domeniu = [int(x) for x in f.readline().split()]
parametrii = [int(x) for x in f.readline().split()]#polinom de grad 2
precizia = int(f.readline()) #nr de zecimale
prob_recombinare = float(f.readline())
prob_mutatie = float(f.readline())
nr_etape = int(f.readline())
etilist = False
tip_mutatie = 1

#intreb daca se tine cont de criteriul etilist sau  nu:
print("Se tine cont de criteriul etilist? Introduceti da sau nu. ")
optiune = str(input())
if optiune.lower() == "da":
    etilist  = True

#intreb ce tip de mutatie e rara sau pentru fiecare gena:
print("Ce tip de mutatie este? Introduceti 1 pentru rara si 2 pt fiecare gena")
optiune = int(input())
if optiune == 2:
    tip_mutatie = 2

#cu ajutorul parametrilor polinomului de grad 2 creez functia explicit de tipul x^2 -x +2 pentru a o afisa
functie = "{}x^2".format(parametrii[0])
if(parametrii[1] >0):
    functie += " +{}x".format(parametrii[1])
elif(parametrii[1] == 0):
    functie += str(parametrii[2])
else:
    functie += " {}x".format(parametrii[1])

if parametrii[1] != 0:
    if parametrii[2] >0:
        functie += " +" + str(parametrii[2])
    elif parametrii[2] <0:
        functie += " " + str(parametrii[2])

#inplementez functia f pt a afisa valorile din cromozomi
def functiepoz(x):
    return parametrii[0] * x**2 + parametrii[1]*x + parametrii[2]


def cautare_binara(arr, threshold):
    low = 0
    high = len(arr)

    while low < high:
        mid = (low + high) // 2


        if arr[mid] == threshold:
            return mid
        elif arr[mid] < threshold and mid != low:
            low = mid
        elif arr[mid] > threshold and mid != high:
            high = mid
        else:
            #se termina cu indexul pointing catre primul element mai mare decat minimul
            high = low = low + 1

    return low


def afiseaza_generatie(generatie):
    for i in range(0, nr_populatie):
        val_intreaga = util.ba2int(generatie[i])
        val_intreaga = ((domeniu[1] - domeniu[0]) / (2 ** lungime_cromozom - 1)) * val_intreaga + domeniu[0]
        g.write("{}: {}, x={}, f={}\n".format(i + 1, generatie[i].to01(), val_intreaga, functiepoz(val_intreaga)))


def calculeaza_max(generatie):
    f_max = 0
    cromozom_max = 0
    val_intreaga = 0
    performanta = 0

    for i in range(0, nr_populatie):
        val_intreaga = util.ba2int(generatie[i])
        val_intreaga = ((domeniu[1] - domeniu[0]) / (2 ** lungime_cromozom - 1)) * val_intreaga + domeniu[0]
        performanta += functiepoz(val_intreaga)
        if f_max < functiepoz(val_intreaga):
            f_max = functiepoz(val_intreaga)
            cromozom_max = generatie[i]

    performanta /= nr_populatie
    return val_intreaga,f_max,cromozom_max,performanta




#afisez datele de intrare
print("Datele de intrare sunt:")
print("-----------------------")
print("functia: ", functie)
print("Domeniul de definitie al functiei: ", domeniu)
print("Dimensiunea populatiei: ", nr_populatie)
print("Parametri pentru functia de maximizat (coeficientii polinomului de grad 2): ", parametrii)
print("Precizia cu care se lucrează (cu care se discretizează intervalul): ", precizia)
print("Probabilitatea de recombinare: ", prob_recombinare)
print("Probabilitatea de mutaţie: ", prob_mutatie)
print("Numărul de etape ale algoritmului: ", nr_etape)
if etilist:
    print("Se tine cont de crieteriul etilist")
else:
    print("NU se tine cont de crieteriul etilist")
print("-----------------------")

#calculez lungimea unui cormozom
lungime_cromozom = round(log((domeniu[1] - domeniu[0])*10**precizia,2))
print("\n lungime cromozom: ", lungime_cromozom)

#prima generatie e random
g.write("Populatia initiala:\n")
cromozomi = [] #lista unde voi tine toata populatia
nr = []#lista cu nr din domeniu date de cromozomi
performanta_totala = 0

#generez x cromozomi unde x este dimensiunea populatiei citita din fisierul de intrare
cromozom_max = 0 #pt selectia etilista se adauga automat cromozomul cel mai fit din etapa curenta
f_max = 0

for i in range(0,nr_populatie):
    #fiecare cromozom este un bitarray => un array cu lungimea calculata antetrior de 0 sau 1 biti selectati random
    cromozom = bitarray()
    for j in range(0,lungime_cromozom):# de la 0 la lungime_cromozom-1
        cromozom.append(randint(0,1))
    cromozomi.append(cromozom)
    val_intreaga = util.ba2int(cromozom) #transform reprezentarea din baza 2 intr-o reprezentare in baza 10 / true => face compelemntul nr binar dat si ajunge la 2**N, un nr mereu pozitiv
    val_intreaga = ((domeniu[1] - domeniu[0])/(2**lungime_cromozom - 1))*val_intreaga + domeniu[0] #daca val_intreaga depaseste domeniul de definitie al functiei il aduc la domeniul de definitie
    performanta_totala += functiepoz(val_intreaga)
    nr.append(val_intreaga)
    g.write("{}: {}, x={}, f={}\n".format(i+1, cromozom.to01(), val_intreaga, functiepoz(val_intreaga)))
    if f_max < functiepoz(val_intreaga):
        cromozom_max = cromozom #aleg cromozomul maxim din populatia initiala
        f_max = functiepoz(val_intreaga)

if etilist:
    print("Cel mai fit cromozom: ", cromozom_max.to01())


g.write("\nProbabilitati de selectie:\n")
probabilitati = []
for i in range(0,nr_populatie):
    g.write("cromozom {} probabilitate {}\n".format(i+1, functiepoz(nr[i])/performanta_totala))
    probabilitati.append(functiepoz(nr[i])/performanta_totala)
#print(cromozomi)


g.write("\nIntervale probabilitati de selectie:\n")

intervale = []
for i in range(0,nr_populatie):
    interval_selectie = 0
    for j in range(0,i):
        interval_selectie += probabilitati[j]#calculez intevalele de selectie qi ca suma din probabilitatile de selectie

    intervale.append(interval_selectie)
    g.write(str(interval_selectie) + " ")

g.write("\n")

selectati = []


for i in range(0,nr_populatie):#ruleta invartita de nr_populatie ori
    #se genereaza variabila aleatoarea uniform pe intervalul [0,1)
    u = uniform(0,1)
    #se cauta binar u in vectorul de intervale de selectie intervale
    #vectorul de intervale e deja sortat deci se poate aplica cautarea binara pe el
    indice = cautare_binara(intervale, u)
    if(indice >= nr_populatie):
        #print("indice: " + str(indice))
        indice -= 1
    selectati.append(cromozomi[indice])

    g.write("u={} selectam cromozonul {}\n".format(u,indice+1))



g.write("\nDupa selectie:\n")#cromozomii care participa la recombinare
afiseaza_generatie(selectati)#voi abtine tot nr_populatie cromozomi (unii se repeta)


g.write("\nProbabilitatea de recombinare: {}\n".format(prob_recombinare))
recombinari = []
nu_participa = []
#se parcug 2 cate 2 cromozomi care vor participa la recombinare(cromozomii selectati)
#daca  e nr de cromozomi impar ii grupez pe ultimii 3
for i in range(0, nr_populatie):#pentru fiecare cromozom din selectie
    #generez un u uniform
    u = uniform(0,1)
    #daca u este mai mic decat prob de recombinare atunci cromozomul participa la recombinare altfel nu
    if(u < prob_recombinare):
        g.write("{}: {} u={}<{} participa\n".format(i+1, selectati[i].to01(), u, prob_recombinare))
        recombinari.append(selectati[i])
    else:
        nu_participa.append(selectati[i])#ii salvez pe cei care nu vor participa la recombinare pentru a afisa noua populatie
        g.write("{}: {} u={}\n".format(i + 1, selectati[i].to01(), u, prob_recombinare))




#incrucisarea sau recombinarea
cromozomi_noi = []
#amestec cromozomii care participa la recombinare
recombinari_initiale = recombinari.copy()#fac o copie pentru a vedea apoi indicii in afisare cromozomul x se combina cu cromozomul y
random.shuffle(recombinari)



#din 2 cromozomi din recombinari se creeaza 2 descendeti noi
#daca avem nr par de cromozomi
g.write("\n")
if len(recombinari) >=2 and len(recombinari)%2 == 0:
    for i in range(0, len(recombinari)-1, 2):#din lista a b c d iar cu pasul doi perechile a,b c,d pe care le recombin
        cromozom1 = recombinari[i]
        cromozom2 = recombinari[i+1]
        g.write("Recombinare dintre cromozomul {} cu cromozomul {}:\n".format(recombinari_initiale.index(cromozom1), recombinari_initiale.index(cromozom2)))
        # generez un punct aleator de rupere
        punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului
        g.write("{} {} punct de rupere {}\n".format(cromozom1.to01(), cromozom2.to01(), punct_rupere))

        # iau prima bucata a cromozomului si a doua bucata si le combin cu bucatile celuilalt cromozom
        b11 = cromozom1[0:punct_rupere]
        b12 = cromozom1[punct_rupere:]
        b21 = cromozom2[0:punct_rupere]
        b22 = cromozom2[punct_rupere:]

        crnou1 = b11 + b22
        cromozomi_noi.append(crnou1)
        crnou2 = b21 + b12
        cromozomi_noi.append(crnou2)
        g.write("Rezultat: {} {}\n".format(crnou1.to01(), crnou2.to01()))

elif len(recombinari) >=2: #daca avem nr impar, diferit de 1 de cromozomi atunci luam ultimii 3 pe care ii recombinam
    for i in range(0, len(recombinari)-3, 2):#ii exclud pe ultimii 3
        cromozom1 = recombinari[i]
        cromozom2 = recombinari[i+1]
        g.write("Recombinare dintre cromozomul {} cu cromozomul {}:\n".format(recombinari_initiale.index(cromozom1), recombinari_initiale.index(cromozom2)))

        # generez un punct aleator de rupere
        punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului
        g.write("{} {} punct de rupere {}\n".format(cromozom1.to01(), cromozom2.to01(), punct_rupere))

        # iau prima bucata a cromozomului si a doua bucata si le combin cu bucatile celuilalt cromozom
        b11 = cromozom1[0:punct_rupere]
        b12 = cromozom1[punct_rupere:]
        b21 = cromozom2[0:punct_rupere]
        b22 = cromozom2[punct_rupere:]

        crnou1 = b11 + b22
        cromozomi_noi.append(crnou1)
        crnou2 = b21 + b12
        cromozomi_noi.append(crnou2)
        g.write("Rezultat: {} {}\n".format(crnou1.to01(), crnou2.to01()))


    #acum ii iau pe ultmii 3 si ii recombin obtinand 3 descendenti
    # generez un punct aleator de rupere

    cromozom1 = recombinari[len(recombinari)-3]
    cromozom2 = recombinari[len(recombinari)-2]
    cromozom3 = recombinari[len(recombinari)-1]
    punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului
    g.write("Recombinare dintre cromozomul {} cu cromozomul {}:\n".format(recombinari_initiale.index(cromozom1), recombinari_initiale.index(cromozom2), recombinari.index(cromozom3)))
    g.write("{} {} {} punct de rupere {}\n".format(cromozom1.to01(), cromozom2.to01(), cromozom3.to01() ,punct_rupere))
    b11 = cromozom1[0:punct_rupere]
    b12 = cromozom1[punct_rupere:]
    b21 = cromozom2[0:punct_rupere]
    b22 = cromozom2[punct_rupere:]
    b31 = cromozom3[0:punct_rupere]
    b32 = cromozom3[punct_rupere:]
    crnou1 = b11 + b22
    crnou2 = b21 + b32
    crnou3 = b31 + b12
    cromozomi_noi.append(crnou1)
    cromozomi_noi.append(crnou2)
    cromozomi_noi.append(crnou3)
    g.write("Rezultat: {} {} {}\n".format(crnou1.to01(), crnou2.to01(), crnou3.to01()))

# populatia dupa recombinare => sunt cromozomii recombinati si cei din populatia initiala care nu au intrat in etapa de recombinare
if len(recombinari) == 1:
    cromozomi_noi = recombinari
    cromozomi_noi.extend(nu_participa)
else:
    cromozomi_noi.extend(nu_participa)


print(len(cromozomi_noi))
g.write("\nDupa recombinare:\n")
afiseaza_generatie(cromozomi_noi)


contor = 1
if tip_mutatie == 2:
    g.write("\nProbabilitate de mutatie pentru fiecare gena {}\n".format(prob_mutatie))
    g.write("\nAu fost modificati cromozomii:\n")
    for cromozom in cromozomi_noi:  # mutatie pe fiecare gena
        for j in range(0, len(cromozom)):
            # se genereaza variabila aleatoarea uniform pe intervalul [0,1)
            u = uniform(0, 1)
            if u < prob_mutatie:  # se schimba o gena j in complementul ei
                g.write("{}\n".format(contor + 1))
                if (cromozom[j] == 0):
                    cromozom = cromozom[0:j] + "1" + cromozom[j + 1:]
                else:
                    cromozom = cromozom[0:j] + "0" + cromozom[j + 1:]
        contor += 1


else:
    g.write("\nProbabilitatea de mutatie rara {}\n".format(prob_mutatie))
    g.write("\nAu fost modificati cromozomii:\n")
    for i in range(0, len(cromozomi_noi)):  # mutatie rara
        # se genereaza variabila aleatoarea uniform pe intervalul [0,1)
        u = uniform(0, 1)
        if u < prob_mutatie:  # se schimba o gena generata aleator in complementul ei
            gena = random.randint(0, len(cromozomi_noi[i]))  # un nr intreg intre 0 si lungimea cromozomului
            g.write("{}\n".format(i + 1))
            if (cromozomi_noi[i] == 0):
                cromozomi_noi[i] = cromozomi_noi[i][0:gena] + "1" + cromozomi_noi[i][gena + 1:]
            else:
                cromozomi_noi[i] = cromozomi_noi[i][0:gena] + "0" + cromozomi_noi[i][gena + 1:]




g.write("\nDupa mutatie:\n")
#daca setine cont de crieteriul etilist, includ cromozomul cel mai fit in generatia urmatoare daca nu a fost deja selectat
if etilist and cromozom_max not in cromozomi_noi:
    cromozomi_noi[0] = cromozom_max

afiseaza_generatie(cromozomi_noi)

g.write("\nEvolutia maximului\n")


'''_Grafic_'''
#graficul pentru evolutia maximului
#axa x e domeniu adica numarul polulatiei
x = []
#axa y e codomeniu adica valoarea maxima atinsa in acea populatie f(val_maxima)
y = []

plt.rcParams['figure.figsize'] = [15,10]


'''_Repeta algoritm_'''
for k in range(0,nr_etape-1):#se repeta algormitmul un nr de etape
    cromozomi = cromozomi_noi
    print(len(cromozomi))
    val_intreaga,f_max, cromozom_max, performanta_totala = calculeaza_max(cromozomi) #calculez maximul pt generatia curenta
    g.write("maximul: {}, performanta totala: {}". format(f_max, performanta_totala))#afisez maximul
    g.write("\n")

    #adaug valorile maxime in axa y pentru grafic
    x.append(k+1)
    y.append(f_max)


    probabilitati = []
    for i in range(0, nr_populatie):
        probabilitati.append(functiepoz(nr[i]) / performanta_totala)


    intervale = []
    for i in range(0, nr_populatie):
        interval_selectie = 0
        for j in range(0, i):
            interval_selectie += probabilitati[
                j]  # calculez intevalele de selectie qi ca suma din probabilitatile de selectie

        intervale.append(interval_selectie)

    selectati = []

    for i in range(0, nr_populatie):  # ruleta invartita de nr_populatie ori
        # se genereaza variabila aleatoarea uniform pe intervalul [0,1)
        u = uniform(0, 1)
        # se cauta binar u in vectorul de intervale de selectie intervale
        # vectorul de intervale e deja sortat deci se poate aplica cautarea binara pe el
        indice = cautare_binara(intervale, u)
        if (indice >= nr_populatie):
            #print("indice: " + str(indice))
            indice -= 1
        selectati.append(cromozomi[indice])

    recombinari = []
    nu_participa = []
    # se parcug 2 cate 2 cromozomi care vor participa la recombinare(cromozomii selectati)
    # daca  e nr de cromozomi impar ii grupez pe ultimii 3
    for i in range(0, nr_populatie):  # pentru fiecare cromozom din selectie
        # generez un u uniform
        u = uniform(0, 1)
        # daca u este mai mic decat prob de recombinare atunci cromozomul participa la recombinare altfel nu
        if (u < prob_recombinare):
            recombinari.append(selectati[i])
        else:
            nu_participa.append(
                selectati[i])  # ii salvez pe cei care nu vor participa la recombinare pentru a afisa noua populatie

    # incrucisarea sau recombinarea
    cromozomi_noi = []
    # amestec cromozomii care participa la recombinare
    recombinari_initiale = recombinari.copy()  # fac o copie pentru a vedea apoi indicii in afisare cromozomul x se combina cu cromozomul y
    random.shuffle(recombinari)

    # din 2 cromozomi din recombinari se creeaza 2 descendeti noi
    # daca avem nr par de cromozomi
    if len(recombinari) >= 2 and len(recombinari) % 2 == 0:
        for i in range(0, len(recombinari) - 1,
                       2):  # din lista a b c d iar cu pasul doi perechile a,b c,d pe care le recombin
            cromozom1 = recombinari[i]
            cromozom2 = recombinari[i + 1]

            # generez un punct aleator de rupere
            punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului

            # iau prima bucata a cromozomului si a doua bucata si le combin cu bucatile celuilalt cromozom
            b11 = cromozom1[0:punct_rupere]
            b12 = cromozom1[punct_rupere:]
            b21 = cromozom2[0:punct_rupere]
            b22 = cromozom2[punct_rupere:]

            crnou1 = b11 + b22
            cromozomi_noi.append(crnou1)
            crnou2 = b21 + b12
            cromozomi_noi.append(crnou2)

    elif len(
            recombinari) >= 2:  # daca avem nr impar, diferit de 1 de cromozomi atunci luam ultimii 3 pe care ii recombinam
        for i in range(0, len(recombinari) - 3, 2):  # ii exclud pe ultimii 3
            cromozom1 = recombinari[i]
            cromozom2 = recombinari[i + 1]

            # generez un punct aleator de rupere
            punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului

            # iau prima bucata a cromozomului si a doua bucata si le combin cu bucatile celuilalt cromozom
            b11 = cromozom1[0:punct_rupere]
            b12 = cromozom1[punct_rupere:]
            b21 = cromozom2[0:punct_rupere]
            b22 = cromozom2[punct_rupere:]

            crnou1 = b11 + b22
            cromozomi_noi.append(crnou1)
            crnou2 = b21 + b12
            cromozomi_noi.append(crnou2)

        # acum ii iau pe ultmii 3 si ii recombin obtinand 3 descendenti
        # generez un punct aleator de rupere

        cromozom1 = recombinari[len(recombinari) - 3]
        cromozom2 = recombinari[len(recombinari) - 2]
        cromozom3 = recombinari[len(recombinari) - 1]
        punct_rupere = random.randint(0, len(cromozom1))  # un nr intreg intre 0 si lungimea cromozomului

        b11 = cromozom1[0:punct_rupere]
        b12 = cromozom1[punct_rupere:]
        b21 = cromozom2[0:punct_rupere]
        b22 = cromozom2[punct_rupere:]
        b31 = cromozom3[0:punct_rupere]
        b32 = cromozom3[punct_rupere:]
        crnou1 = b11 + b22
        crnou2 = b21 + b32
        crnou3 = b31 + b12
        cromozomi_noi.append(crnou1)
        cromozomi_noi.append(crnou2)
        cromozomi_noi.append(crnou3)

    # populatia dupa recombinare => sunt cromozomii recombinati si cei din populatia initiala care nu au intrat in etapa de recombinare
    if len(recombinari) == 1:
        cromozomi_noi = recombinari
        cromozomi_noi.extend(nu_participa)
    else:
        cromozomi_noi.extend(nu_participa)

    contor = 1
    if tip_mutatie == 2:
        for cromozom in cromozomi_noi:  # mutatie pe fiecare gena
            for j in range(0, len(cromozom)):
                # se genereaza variabila aleatoarea uniform pe intervalul [0,1)
                u = uniform(0, 1)
                if u < prob_mutatie:  # se schimba o gena j in complementul ei
                    if (cromozom[j] == 0):
                        cromozom = cromozom[0:j] + "1" + cromozom[j + 1:]
                    else:
                        cromozom = cromozom[0:j] + "0" + cromozom[j + 1:]
            contor += 1


    else:
        for i in range(0, len(cromozomi_noi)):  # mutatie rara
            # se genereaza variabila aleatoarea uniform pe intervalul [0,1)
            u = uniform(0, 1)
            if u < prob_mutatie:  # se schimba o gena generata aleator in complementul ei
                gena = random.randint(0, len(cromozomi_noi[i]))  # un nr intreg intre 0 si lungimea cromozomului
                if (cromozomi_noi[i] == 0):
                    cromozomi_noi[i] = cromozomi_noi[i][0:gena] + "1" + cromozomi_noi[i][gena + 1:]
                else:
                    cromozomi_noi[i] = cromozomi_noi[i][0:gena] + "0" + cromozomi_noi[i][gena + 1:]

    print(k, len(recombinari), len(nu_participa),len(cromozomi_noi))
    # daca setine cont de crieteriul etilist, includ cromozomul cel mai fit in generatia urmatoare daca nu a fost deja selectat
    if etilist and cromozom_max not in cromozomi_noi:
        cromozomi_noi[0] = cromozom_max













#tipul graficului
plt.plot(x, y, color='green', linestyle='dashed', linewidth=3,
         marker='o', markerfacecolor='blue', markersize=10)


#dau deumiri axelor
plt.xlabel("axa x")
# naming the y axis
plt.ylabel("axa y")

#titlu pt grafic
plt.title("Evolutia maximului:")



plt.show()

