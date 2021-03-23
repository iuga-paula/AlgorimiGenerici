#algoritm genetic pentru determinarea maximului unei funcţii pozitive pe un domeniu dat
from math import log
from bitarray import bitarray, util
from random import randint, choice, uniform

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

        #print("low = ", low, " mid = ", mid, " high = ", high)

        if arr[mid] == threshold:
            return mid
        elif arr[mid] < threshold and mid != low:
            low = mid
        elif arr[mid] > threshold and mid != high:
            high = mid
        else:
            # terminate with index pointing to the first element greater than low
            high = low = low + 1

    return low

#def recombinare(cromozom1, cromozom2):

#afisez datele de intrare
print("Datele de intrare sunt:")
print("-----------------------")
print("functia: ", functie)
print("Domeniul de definitie al functiei: ", domeniu)
print("Dimensiunea populatiei: ", nr_populatie)
print("Parametri pentru functia de maximizat (coeficientii polinomului de grad 2): ", parametrii)
print("precizia cu care se lucrează (cu care se discretizează intervalul): ", precizia)
print("probabilitatea de recombinare: ", prob_recombinare)
print("probabilitatea de mutaţie: ", prob_mutatie)
print("numărul de etape ale algoritmului: ", nr_etape)
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
        cromozom_max = cromozom.to01()
        f_max = functiepoz(val_intreaga)

print("Cel mai fit cromozom: ", cromozom_max)

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
        interval_selectie += probabilitati[j]#calculez intevalele de lectie qi ca suma din probabilitatile de selectie

    intervale.append(interval_selectie)
    g.write(str(interval_selectie) + " ")

g.write("\n")

selectati = []


for i in range(0,nr_populatie):
    #se genereaza variabila aleatoarea uniform pe intervalul [0,1)
    u = uniform(0,1)
    #se cauta binar u in vectorul de intervale de selectie intervale
    #vectorul de intervale e deja sortat deci se poate aplica cautarea binara pe el
    indice = cautare_binara(intervale, u)
    if(indice >= nr_populatie):
        indice -= 1
    selectati.append(cromozomi[indice])

    g.write("u={} selectam cromozonul {}\n".format(u,indice+1))



g.write("\nDupa selectie:\n")#cromozomii care participa la recombinare
for i in range(0,nr_populatie):
    val_intreaga = util.ba2int(selectati[i])
    val_intreaga = ((domeniu[1] - domeniu[0]) / (2 ** lungime_cromozom - 1)) * val_intreaga + domeniu[0]
    g.write("{}: {}, x={}, f={}\n".format(i + 1, selectati[i].to01(), val_intreaga, functiepoz(val_intreaga)))

g.write("\nProbabilitatea de recombinare: {}\n".format(prob_recombinare))
recombinari = []
#se parcug 2 cate 2 cromozomi care vor participa la recombinare(cromozomii selectati)
#daca  e nr de cromozomi impar il ignor pe ultimul
for i in range(0, nr_populatie):#pentru fiecare cromozom din selectie
    #generez un u uniform
    u = uniform(0,1)
    #daca u este mai mic decat prob de recombinare atunci cromozomul participa la recombinare altfel nu
    if(u < prob_recombinare):
        g.write("{}: {} u={}<{} participa\n".format(i+1, selectati[i].to01(), u, prob_recombinare))
        recombinari.append(selectati[i])
    else:
        g.write("{}: {} u={}\n".format(i + 1, selectati[i].to01(), u, prob_recombinare))






