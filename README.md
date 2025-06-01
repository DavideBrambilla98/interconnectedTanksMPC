# Progettazione di un MPC per il controllo del livello dell’acqua in quattro serbatoi interconnessi

Il quadruple tank process è un sistema multivariabile composto da serbatoi interconnessi che viene spesso utilizzato come benchmark per controllori MPC, date la non-linearità della sua dinamica e la necessità di definire dei vincoli sulle sue variabili di stato.

## 1. Descrizione del problema
Si controlli il livello dell’acqua nei quattro serbatoi interconnessi riportati in figura 1.

![image](https://github.com/user-attachments/assets/733e57e7-1867-48fb-be24-cd1afd9983ff)

Figure 1: Schema dei quattro serbatoi interconnessi.

Si ipotizzi che le due pompe del sistema siano controllate tramite le tensioni v1 e v2 e che la configurazione
delle valvole responsabili della divisione dei flussi d’acqua tra i serbatoi sia lasciata inalterata
durante il funzionamento dell’impianto.
La dinamica del sistema è definita dalle equazioni:

<img width="169" alt="image" src="https://github.com/user-attachments/assets/45160abd-d3cb-47e5-a92f-5647410b7440" />

dove:

• hi [cm] è il livello dell’acqua nell’i-esimo serbatoio;

• Ai [cm2] è la sezione dell’i-esimo serbatoio;

• ai [cm2] è la sezione del foro presente nell’i-esimo serbatoio;

• g [cm/s2] è l’accelerazione di gravità;

• kivi [cm3/s] è il flusso d’acqua generato dall’i-esima pompa, controllata tramite la tensione vi [V];

• γi [/] definisce come il flusso generato dall’i-esima pompa viene diviso.

I valori dei parametri presenti nelle precedenti equazioni sono riportati nella seguente tabella:

<img width="170" alt="image" src="https://github.com/user-attachments/assets/7297dee6-2675-4708-9c8a-f44a5caa6bcc" />

Le tensioni vi applicate alle due pompe sono non-negative e non possono essere superiori ai 4.5V.
I livelli dell’acqua nei serbatoi devono sempre rimanere nel range [0.5cm, 20cm].
Gli ingressi controllati del sistema sono le tensioni applicate alle due pompe.


## 2. Obiettivi
L’obiettivo del lavoro è progettare un controllore MPC in grado di portare il sistema dalla condizione
iniziale (h1, h2, h3, h4) = (1.3767, 2.2772, 0.8386, 0.5604) all’equilibrio (7.8253, 18.7323, 3.3545, 7.8801),
rispettando sempre i vincoli, e simulare il funzionamento del sistema in anello chiuso.

Si confrontino le prestazioni del controllore progettato al variare di:

• ingredienti terminali (vincolo terminale di uguaglianza vs costo terminale e vincolo terminale di
disuguaglianza);

• Q ed R del costo quadratico;

• orizzonte di predizione;

• tempo di campionamento ts ≥ 1s.
