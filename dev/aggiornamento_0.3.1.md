# Aggiornamento di fewclust alla versione 0.3.1

Preparato il 2026-11-02 nella sessione parallela sul pacchetto (ramo `package`, worktree `~/MMM-package`), dopo la decisione del PI sugli esiti costanti (decisione D-82 del ramo del manoscritto, `docs/decisions/D-82.md` su `morris-v4`: "Un dataset con esito costante vale nessun risultato per ogni metodo"). Il ramo del manoscritto (`morris-v4`, `~/MMM`) si legge soltanto. Stessa forma dell'aggiornamento precedente (`dev/aggiornamento_0.3.0.md`).

**Stato di questo file.** Le sezioni da "Obiettivo" a "Test esistenti che cambiano" sono scritte e committate *prima* di toccare il codice del pacchetto, insieme agli script dei controlli in `dev/checks_0.3.1/` e ai test di `tests/testthat/` riscritti sulle nuove attese. La sezione "Esiti" si aggiunge alla fine, con i numeri letti dai file in `dev/checks_0.3.1/out/`.

## Obiettivo

Portare fewclust dalla 0.3.0 (commit `22e618e`) a una 0.3.1 la cui scelta automatica sia quella dello studio corretto. Nella 0.3.0 la scelta era la GEE con correlazione di lavoro AR(1), la correzione di Fay e Graubard e il riferimento t. Lo studio contava però come risultato usabile della GEE 9.646 dataset con esito tutto 0 o tutto 1. Tolti quelli, la GEE AR(1) con Fay e Graubard e t è chiaramente sopra il 7,5% di errore di I tipo in 2 punti di disegno (scenario M2, 10 cluster, 2 occasioni, meccanismi GLMM-native e GLMM-AR1). La regola dello studio (strategie con un intervallo; si escludono quelle chiaramente sopra 7,5% in qualche punto; fra le altre si sceglie la quota di risultati usabili più alta) sceglie allora il modello marginalizzato multilivello con il controllo di convergenza e la t: 0 punti chiaramente sopra, quota usabile 90,9%. I calcoli che il paper ha preso dal pacchetto (`gee_vcov()`, motore del modello marginalizzato) non cambiano.

## Da dove vengono i numeri

Si leggono dal ramo del manoscritto, mai ricopiati a mano. I file sono nella copia di lavoro di `~/MMM` (step S9 aperto, `pi_approval: pending`): **non sono ancora committati**. Il pacchetto si riferisce alla versione con questi md5, verificati all'inizio e alla fine dei controlli:

| file | contenuto | md5 |
|---|---|---|
| `tables/V5_S7f_table4.csv` | tabella per strategia (Tabella 5 del paper v6; nel pacchetto `table4`) | `bae9d66a6dda8790e8cdd1537a062888` |
| `tables/V5_S7f_h2h_boot.csv` | confronto diretto fra le due strategie per regione (Figura 1a), classi | `14eeef1dbc4c46cd9692b2420331d3bf` |
| `tables/V5_S7f_rule_dp.csv` | per punto di disegno e strategia: chiaramente sopra / sotto, quota usabile | `494a22502e437a372f7210a7bb3f39cf` |
| `tables/V5_S7f_cp_S7d_dp_rates.csv` | per punto di disegno: dataset usabili, rifiuti sotto H0, tasso | `e3b6711c6f423493081ce98f422444d9` |

Se uno di questi md5 cambia prima della chiusura, i controlli falliscono e le costanti si rigenerano.

## Cosa cambia (piano)

1. **Scelta automatica.** `auto_cluster_fit()` e `cluster_strategies()` scelgono il modello marginalizzato multilivello con il controllo di convergenza e la t su N - p gradi di libertà (`"mmm"`), a ogni numero di cluster e occasioni. L'alternativa è la GEE AR(1) con Fay e Graubard e t (`"gee_ar1_fg"`): il piano la nomina in `$alternative` e `auto_cluster_fit()` la stima, registrando la ricaduta nel percorso, se il modello marginalizzato non restituisce un risultato usabile. Con 2-3 occasioni e al massimo 25 cluster `flip_test()` resta il complemento (`$complement = "flip_test"`), come nella 0.3.0. Se l'utente cambia la correlazione di lavoro in `fewclust_control(corstr = )`, l'alternativa la usa (per esempio `"gee_exch_fg"`) e una cautela dice che la regola dello studio usa AR(1). Se falliscono entrambe, il rifiuto suggerisce `flip_test()`. La cautela sui buchi nelle serie (un soggetto con un'occasione mancante in mezzo) vale ora anche per il modello marginalizzato, il cui termine di transizione lega le osservazioni consecutive del soggetto: nella 0.3.0 valeva solo per la correlazione AR(1), che era la scelta automatica.
2. **Cautela sulla GEE AR(1) con Fay e Graubard a 10 cluster e 2 occasioni.** Quando l'utente sceglie `"gee_ar1_fg"`, o la ricaduta la usa, con meno di 15 cluster e al più 2 occasioni, una cautela cita i 2 punti di disegno dello studio fuori taratura: meccanismo, scenario, rifiuti sotto H0 su dataset usabili, tasso e quota usabile, letti da `V5_S7f_rule_dp.csv` e `V5_S7f_cp_S7d_dp_rates.csv`. A ogni dimensione, come per le altre strategie, la cautela generale cita il conteggio della tabella (2 punti su 2.376).
3. **Note per regione dal nuovo confronto diretto.** Le note della 0.3.0 basate sul vecchio confronto ("con 2 occasioni e 15 o più cluster il modello marginalizzato era un po' più vicino al 5%") si tolgono. Al loro posto, per la cella della mappa della Figura 1a in cui cade il disegno (cluster: 10, 15-25, oltre 25; occasioni: 2, 3, 4 o più), se la classe di `mad_gee_minus_mmm` in `V5_S7f_h2h_boot.csv` dice che una strategia è più vicina al 5% dell'altra oltre il margine (classe `GEE` o `MMM`), una nota lo dice con la distanza media dal 5% delle due strategie e il numero di punti. Con il file attuale questo accade solo a 10 cluster con 3 e con 4 o più occasioni (a favore della GEE); le classi `similar` e `about the margin` non danno nota. La nota compare sotto la scelta automatica e quando l'utente sceglie `"mmm"` o `"gee_ar1_fg"`.
4. **Esito costante.** Nella 0.3.0 (verificato con la 0.3.0 installata, prima di questo piano) `auto_cluster_fit()` e `cluster_strategies()` rifiutano già un dataset con esito costante con la regola "D2: constant outcome", per ogni strategia, anche se richiesta dall'utente; `compare_strategies()` parte da un fit, quindi non arriva mai a un tale dataset. `flip_test()` invece restituisce `NA` con il motivo "score numerically zero", che non dice perché. Nella 0.3.1 `flip_test()` controlla l'esito prima del calcolo (dopo aver estratto i segni, così il flusso casuale non cambia) e restituisce `NA` con il motivo "constant outcome (all 0 or all 1): no usable result, as for every strategy". `gee_vcov()` non è una strategia (calcola la varianza di un fit `geeglm` fatto dall'utente) e non cambia, per il controllo di invarianza.
5. **Costanti rigenerate.** `R/constants.R` si scrive con `dev/make_constants.R`, che legge i quattro file sopra: la tabella per strategia (tutte le 19 righe), la strategia scelta dalla regola (colonna `chosen`), i punti chiaramente sopra 7,5% della GEE AR(1) con Fay e Graubard e la mappa della Figura 1a. Nessun numero si scrive a mano. Il testo per l'utente chiama la tabella "Table 5 of the study" (numerazione del paper v6); nel codice resta `table4`.
6. **Documentazione.** `DESCRIPTION` versione 0.3.1 e descrizione riscritta; voce 0.3.1 in `NEWS.md`; `README.md`, le due vignette, `inst/CITATION` (manoscritto v6) e la documentazione delle funzioni (`auto_cluster_fit`, `cluster_strategies`, `flip_test`, `fewclust_control`, pacchetto) riscritte sulla nuova raccomandazione; `man/` rigenerato con roxygen2.

### Una scelta per il PI ancora aperta

La condizione di coerenza del controllo di convergenza (`fewclust_control(check_coherence = TRUE)`, un fit con il trend la cui log-verosimiglianza è sotto quella del fit senza trend si tratta come non convergente) **resta attiva per default, come nella 0.3.0**. Non è nella regola del paper: nella 0.3.0 cambiava il fit riportato in 13 dataset su 796 dello scenario con deviazione standard latente crescente. Se tenerla è una decisione del PI ancora aperta dalla 0.3.0, e questo aggiornamento non la prende. Con la 0.3.1 però pesa di più: il modello marginalizzato diventa la scelta automatica, quindi il default del pacchetto può differire dalla regola esatta del paper proprio nella strategia raccomandata. Il controllo 3 sotto confronta il pacchetto con lo studio usando la regola del paper (`check_coherence = FALSE`) e riporta, come descrizione, in quanti dataset il default dà un fit diverso.

### Scelte tecniche dell'agente

- La cautela del punto 2 si attiva sotto i 15 cluster (lo studio ha simulato 10 e poi 15) e con al più 2 occasioni. Fra 11 e 14 cluster lo studio non ha dati: la cautela vale come prudenza.
- Le bande della mappa (punto 3) seguono le regioni del file (`map|10|...`, `map|15-25|...`, `map|>25|...`); un disegno con 11-14 cluster cade nella banda "10", con 26-34 nella banda "oltre 25".
- L'avviso del modello marginalizzato a 10 cluster ("conservative ... 18 of the 156 design points with 10 clusters") resta, senza l'inciso "(the GEE strategy was closer to 5% there)": con 2 occasioni a 10 cluster il confronto è ora "about the margin"; con 3 o più ci pensa la nota del punto 3.

## Controlli (cancelli)

Ognuno è binario; se uno fallisce, la 0.3.1 non si chiude. Script in `dev/checks_0.3.1/` (adattati da `dev/checks_0.3.0/`), lanciati da `run_checks.sh`; risultati in `dev/checks_0.3.1/out/`. Libreria di sviluppo `~/R/fewclust_dev_lib`; `~/R/s7d_lib` non si tocca. Al massimo 4 core. I file temporanei (la 0.3.0 installata per il confronto, i file del controllo di invarianza, la cartella di `R CMD check`) stanno in `~/MMM-package/.tmp031/`, non in `/tmp`, e si cancellano alla fine.

1. **Test a inversione di segno riprodotto** (`c1_study_repro.R`, stesso insieme di 126 dataset della 0.3.0). Come nella 0.3.0: p-value uguale a quello salvato entro 1e-12 oppure entrambi non calcolabili con lo stesso motivo; in più, sui dataset con esito costante il pacchetto deve dire "constant outcome" e il p-value salvato deve essere `NA` (qualunque il motivo salvato). Passa se coincidono tutti, con almeno 20 dataset e dati identici al run salvato.
2. **L'alternativa GEE riprodotta** (stesso script). Sugli stessi dataset, chiedendo esplicitamente `strategy = "gee_ar1_fg"`, stima ed errore standard di `time:trt` uguali a quelli salvati (`ar1_est`, `ar1_se_fg`) con tolleranza relativa 1e-8, dove il pacchetto non rifiuta e lo studio ha un errore standard finito. Passa se coincidono tutti, su almeno 20 confronti. Lo script registra anche la strategia scelta da `auto_cluster_fit()`.
3. **La strategia raccomandata riprodotta** (`c9_mmm_repro.R`, nuovo). Dataset rigenerati dai semi con il codice di `morris-v4`: scenario M3 (deviazione standard latente con trend nel tempo, l'unico scenario in cui lo studio ha salvato il modello marginalizzato con il controllo per dataset con lo stesso modello del pacchetto: media `binY ~ time * trt`, deviazione standard latente `~ time`, transizione `~ 1`), meccanismi GenBinaryY, GLMM-native, GLMM-AR1, N in {10, 25, 65}, T in {2, 4, 8}, sotto H0 e H1, replicati 1-3 (162 dataset). Riferimento: `cc_est`, `cc_se` e `used` di `results_v4/S6cc/rows/` (MMLB::mm, q = 20, controllo di convergenza del paper). Su ciascuno `auto_cluster_fit(binY ~ time * trt, group = id, time = time, control = fewclust_control(check_coherence = FALSE))` deve scegliere `"mmm"`, seguire lo stesso ramo del controllo (`used` = B se e solo se "constant fit used") e dare stima ed errore standard di `time:trt` uguali con `all.equal` a tolleranza 1e-6 (stima) e 1e-5 (errore standard), le tolleranze di `test-mmm-check.R`. Non si confrontano numericamente, e si contano: i dataset non generati o con disegno degenere, quelli che il pacchetto rifiuta (per esempio esito costante), quelli dove lo studio non ha un fit con il termine, e quelli dove il fit riportato dallo studio ha codice `nlm` diverso da 1 e 2 (il punto d'arresto dipende dall'implementazione, correzione 1 della 0.3.0; lì si confronta solo il ramo). Passa se coincidono tutti i confronti, almeno 20, e se i dati rigenerati hanno lo stesso numero di righe e di eventi dei run salvati (`n_obs`, `sum_y` di `results_v4/S7d2/rows/`). Descrizione, non cancello: in quanti dataset il default del pacchetto (`check_coherence = TRUE`) riporta un fit diverso.
4. **Regola sulla griglia 3 × 3** (`c2_rule_fig1.R` e `tests/testthat/test-decision.R`). Per N in {10, 25, 145} e T in {2, 3, 10}, `cluster_strategies()` deve dare: strategia uguale alla riga `chosen` di `V5_S7f_table4.csv` (`mmm_cc`, cioè `"mmm"`), riferimento t, alternativa `"gee_ar1_fg"`, complemento `"flip_test"` se e solo se T ≤ 3 e N ≤ 25, nota di regione se e solo se la classe della cella della mappa in `V5_S7f_h2h_boot.csv` è `GEE` o `MMM`, nessuna nota "slightly closer to 5%". Passa se tutte e 9 le celle coincidono.
5. **`R CMD check --as-cran`** (`c7_rcmdcheck.sh`): nessun errore e nessun warning; ogni NOTE motivato in `cran-comments.md`.
6. **Test del pacchetto** (`testthat`): tutti passano. I test delle correzioni (`test-corrections.R`), del motore del modello marginalizzato (`test-mmm.R`) e del controllo di convergenza (`test-mmm-check.R`) restano con le stesse attese numeriche.
7. **I calcoli presi dal paper non cambiano** (`c3_invariance_*.R`). `gee_vcov()` su 100 dataset simulati, tre correlazioni di lavoro, `scale.fix` vero e falso, e il motore del modello marginalizzato su 40 dataset: identici bit per bit fra la 0.3.0 (commit `22e618e`, installata in `~/MMM-package/.tmp031/lib030`) e la 0.3.1.
8. **Esito costante: nessun risultato** (`c8_constant_outcome.R`, nuovo). Dataset costruiti con `simulate_clustered_binary()` per N in {10, 25} e T in {2, 4}, con e senza una covariata, con l'esito messo tutto a 0 e tutto a 1 (16 dataset). Su ciascuno: `cluster_strategies()` rifiuta con la regola "D2: constant outcome"; `auto_cluster_fit()` restituisce un `fewclust_refusal` con quella regola per la scelta automatica e per ognuno dei 19 nomi di strategia accettati (`fewclust:::fc_strategy_names()`: i 7 nomi base e le 12 GEE con la correlazione nel nome); `flip_test()` restituisce p-value `NA`, `computable = FALSE`, motivo che inizia con "constant outcome". Controprova: con un solo esito cambiato il motivo "constant outcome" non compare. Passa se tutto vale su tutti i 16 dataset.
9. **Costanti uguali ai file** (`c5_constants.R`). Ogni valore numerico della tabella per strategia (19 strategie, tutte le colonne numeriche copiate), la strategia scelta, i 2 punti fuori taratura della GEE AR(1) con Fay e Graubard (con dataset usabili, rifiuti, tasso, quota usabile) e le righe della mappa sono uguali ai file (tolleranza relativa 1e-9); gli md5 dei quattro file sono quelli di questo piano; nessun testo del pacchetto cita ancora le frasi della vecchia raccomandazione (`slightly closer to 5%`, `1,896`, "recommended strategy ... is generalized estimating equations"); `DESCRIPTION` ha la versione 0.3.1 e `NEWS.md` si apre con la 0.3.1.
10. **Condizione per citare il pacchetto, adattata** (`c4_paper_condition.R`). Lo stesso codice della 0.3.0 sugli stessi due dataset, con la condizione aggiornata alla nuova regola: la strategia scelta è `"mmm"`. Lo script del ramo del manoscritto (`R/s7f/s7f_fewclust_condition.R`) cerca ancora "ar1" e "fg": va aggiornato nella chat del manoscritto, non qui.

Integrità: `run_checks.sh` verifica gli md5 dei quattro file d'ingresso all'inizio e alla fine (cancello 9). Commit e stato di `~/MMM` si registrano nel log ma non sono un cancello: un'altra sessione lavora su quel ramo e la sua copia di lavoro può cambiare durante i controlli senza che questa sessione la tocchi. `c1` e `c9` confrontano l'md5 dei file di risultati letti prima e dopo.

Il controllo P1 della 0.3.0 (peso della condizione di coerenza, `c6_coherence.R`) non si ripete: il codice del modello marginalizzato non cambia (cancello 7). La descrizione del cancello 3 dà il dato aggiornato sulla nuova scelta automatica.

## Test esistenti che cambiano

- `test-decision.R`: scelta automatica `"mmm"` con alternativa `"gee_ar1_fg"`; la griglia 3 × 3 con le note di regione; conteggi della nuova tabella (per esempio `gee_ar1_kc` 138, `gee_ar1_md` 20, `gee_ind_fg` 5, `gee_exch_fg` 15, `glmm_wald` 515, `glmm_lrt` 412; `"gee_fg"` ora con "2 of 2,376", `"mmm"` senza cautele a 30 cluster); la cautela a 10 cluster e 2 occasioni; la cautela dei disegni brevi con "39% (marginalized model) and 44% (GEE alternative)"; la correlazione di lavoro diversa da AR(1); la cautela del modello marginalizzato a 10 cluster ora anche sotto la scelta automatica.
- `test-fit.R`: la scelta automatica è il modello marginalizzato; la ricaduta va dal modello marginalizzato alla GEE AR(1) con Fay e Graubard, con la sua cautela; il rifiuto quando falliscono entrambe; la GEE con Fay e Graubard si chiede esplicitamente.
- `test-methods.R`, `test-compare.R`, `test-gee.R`, `test-plot.R`: la strategia scelta è `"mmm"`, la GEE è un override.
- `test-flip.R`: esito costante, motivo "constant outcome", senza cambiare il flusso casuale.
- `test-mmm.R`, `test-mmm-check.R`, `test-corrections.R`: nessuna attesa numerica cambia.

## Esiti

Scritti il 2026-11-02 dopo l'ultimo giro completo dei controlli sul codice finale (`run_checks.sh all`, log in `out/run_checks.log`), più la diagnosi del controllo 3 (`c9b_diagnose.R`, scritta dopo il fallimento, non è un cancello). Ogni numero è letto dai file di `dev/checks_0.3.1/out/`.

### In breve: il pacchetto ora sceglie quello che sceglie lo studio corretto?

**A cosa serviva.** Nello studio, i dataset in cui l'esito è tutto 0 o tutto 1 contavano come risultato usabile per la GEE. Tolti quelli, la GEE con correlazione AR(1) e correzione di Fay e Graubard va fuori taratura in due situazioni con 10 soggetti e 2 visite, e la regola dello studio sceglie il modello marginalizzato multilivello con il controllo di convergenza. Il pacchetto sceglieva ancora la GEE: andava allineato, senza cambiare i calcoli che il paper ha preso dal pacchetto.

**Cosa abbiamo trovato.**

1. **Il pacchetto ora sceglie il modello marginalizzato con il controllo di convergenza, e usa la GEE AR(1) con Fay e Graubard solo se il primo non dà un risultato.**
   - Su 9 combinazioni di soggetti e visite la scelta, l'alternativa e il test a inversione di segno come complemento sono quelli dello studio corretto in tutte e 9.
   - Se si usa la GEE con meno di 15 soggetti e 2 visite, un avviso cita le due situazioni dello studio in cui sbagliava troppo: fra 17 e 20 falsi positivi su 100 test, su appena il 6% di dataset con un risultato.
   - Con 10 soggetti e 3 o più visite una nota dice che lì la GEE era più vicina al 5%; la vecchia nota "con 2 visite il modello marginalizzato era un po' migliore" non c'è più, perché nel confronto corretto quella differenza è al limite della soglia di importanza.
2. **Il modello marginalizzato del pacchetto riproduce quello dello studio, salvo 2 dataset su 119 con un risultato usabile: il controllo preregistrato quindi non passa.**
   - Su 162 dataset dello studio rigenerati, 127 sono confrontabili; in 119 il pacchetto dà la stessa stima e lo stesso errore standard.
   - 6 delle 8 differenze sono dataset con pochissimi eventi, dove la stima esplode e né lo studio né il pacchetto hanno un risultato usabile: lì cambia solo un errore standard di migliaia.
   - Le 2 differenze vere: in un dataset il fit di riserva dello studio aveva una stima finita per il termine saggiato ma non per tutti i coefficienti, e il pacchetto, che li chiede tutti finiti, non lo usa; nell'altro l'ottimizzatore si ferma, in entrambe le implementazioni, su un falso massimo, in punti un po' diversi (errore standard 0,12 contro 0,17, mentre il fit senza trend dà 0,31).
3. **La GEE alternativa e il test a inversione di segno riproducono lo studio fino all'ultima cifra.** Stima ed errore standard della GEE coincidono in 51 dataset su 51; il p-value del test coincide in 126 su 126, compresi 9 dataset con esito costante, dove ora il pacchetto dice esplicitamente "esito costante" invece del motivo tecnico che dava prima.
4. **Un dataset con esito costante non dà risultato con nessun metodo del pacchetto.** Su 16 dataset costruiti apposta, tutte le 20 scelte possibili (automatica e 19 strategie) vengono rifiutate, e il test a inversione di segno non dà p-value.
5. **I numeri del paper presi dal pacchetto non cambiano.** Le correzioni del sandwich sono identiche bit per bit alla versione precedente su 600 combinazioni, il modello marginalizzato su 40 dataset; tutti i 344 numeri citati dal pacchetto sono uguali alle tabelle dello studio; il controllo formale per CRAN non dà né errori né avvisi.

**Cosa devi decidere.**

1. **Accettiamo le 2 differenze del modello marginalizzato come spiegate e chiudiamo la versione? Consiglio di sì, con una frase nella documentazione.**
   - Ragioni: riguardano 2 dataset su 119 con un risultato, entrambi in casi in cui l'ottimizzatore non converge davvero; il pacchetto in un caso non dà risultato (prudente), nell'altro riporta, come lo studio, un fit fermo su un falso massimo.
   - Alternativa: allineare il fit di riserva alla regola dello studio (basta che il termine saggiato sia finito). Il pacchetto diventerebbe uguale allo studio in quel caso, ma potrebbe riportare un risultato con altri coefficienti non stimati; serve un altro giro di controlli.
2. **Teniamo attiva per default la condizione in più del controllo di convergenza (quella sulla verosimiglianza)? Ora conta di più, perché il modello marginalizzato è la scelta automatica. Consiglio di tenerla.**
   - Su questi 127 dataset cambia il risultato in 4; in tutti e 4 il fit scartato aveva un errore standard fra un terzo e la metà di quello del fit senza trend. Il secondo caso del punto precedente è uno di questi: la condizione in più lo intercetta.
   - Se la spegniamo, il pacchetto fa esattamente la regola dello studio, compresi questi casi.

**Comunque.** Nella chat del manoscritto bisognerà aggiornare il controllo con cui il paper decide se citare il pacchetto: cerca ancora la GEE e va riscritto per il modello marginalizzato. Le tabelle dello studio da cui vengono i numeri del pacchetto non sono ancora committate nel ramo del manoscritto: se cambiano prima della chiusura di quella sessione, i numeri del pacchetto si rigenerano con lo stesso script.

Appena decidi questi 2 punti, la 0.3.1 è chiusa.

### Cancelli

| Controllo | Criterio | Osservato | Esito |
|---|---|---|---|
| 1. Test a inversione di segno riprodotto | tutti i p-value entro 1e-12 o non calcolabili con lo stesso motivo; esito costante: motivo "constant outcome" e p salvato `NA`; ≥ 20 dataset; dati identici | 126 su 126 (73 con tutte le 2^N assegnazioni, 51 con 1.999 casuali, 49 non calcolabili nel run salvato, di cui 9 con esito costante; 2 non generabili qui e nel run salvato); differenza massima 0; dati identici 126 su 126; md5 invariato | passa |
| 2. L'alternativa GEE riprodotta | stima e SE con tolleranza relativa 1e-8; ≥ 20 confronti | 51 su 51, differenza massima 0 | passa |
| 3. Strategia raccomandata riprodotta | stessa strategia, stesso ramo, stima `all.equal` 1e-6, SE 1e-5; ≥ 20 confronti; dati identici | 119 su 127 confronti (114 numerici); dati identici 162 su 162; md5 invariato | **non passa** |
| 4. Regola sulla griglia 3 × 3 | 9 celle su 9 | 9 su 9 | passa |
| 5. `R CMD check --as-cran` | 0 errori, 0 warning | `Status: 1 NOTE` (validazione HTML del manuale saltata: `tidy` non installato, come nella 0.3.0) | passa |
| 6. Test del pacchetto | tutti passano | 15 file, 706 attese, 0 fallite, 0 errori | passa |
| 7. Calcoli presi dal paper invariati | identici bit per bit fra 0.3.0 e 0.3.1 | `gee_vcov()` 600 insiemi su 600 (64 con messaggio d'errore, identico anch'esso); motore del modello marginalizzato 40 su 40 | passa |
| 8. Esito costante | 16 dataset, rifiuto D2 per 20 scelte, `flip_test()` non calcolabile con motivo "constant outcome"; controprova | 16 su 16; 320 fit su 320 rifiutati con D2; controprova 16 su 16 | passa |
| 9. Costanti uguali ai file | ogni valore uguale; md5 come nel piano; nessun testo vecchio; versione e NEWS | 344 valori su 344 (285 della tabella per strategia, 14 dei 2 punti fuori taratura, 45 della mappa); strategia scelta `mmm_cc`; md5 uguali; 0 righe vecchie; versione 0.3.1 | passa |
| 10. Condizione per citare il pacchetto, adattata | scelta `"mmm"` sui due dataset | 2 su 2 | passa |

Fonti: `out/c1_summary.txt` (1 e 2), `out/c9_summary.txt` e `out/c9b_summary.txt` (3), `out/c2_summary.txt`, `out/c7_summary.txt` e `out/c7_check.log`, `out/tests.log`, `out/c3_summary.txt`, `out/c8_summary.txt`, `out/c5_summary.txt`, `out/c4_summary.txt`. Integrità: gli md5 dei quattro file d'ingresso sono gli stessi all'inizio e alla fine (`out/run_checks.log`: "inputs unchanged during the checks: yes"); il commit di `~/MMM` è rimasto `6b2159e`, i file modificati nella sua copia di lavoro sono passati da 143 a 144 per il lavoro dell'altra sessione.

### Controllo 3: perché non passa

Diagnosi in `out/c9b_summary.txt` e `out/c9b_mismatches.csv`. Dei 162 dataset, 35 non sono confrontati (14 rifiutati dal pacchetto perché l'esito varia in un solo soggetto, 17 per separazione, uno dei quali anche senza fit usabile nello studio, 4 in cui né il modello marginalizzato né la GEE del pacchetto danno un risultato). Nessuno ha esito costante. Delle 8 differenze sui 127 confronti:

- **6 non hanno un risultato usabile né nello studio né nel pacchetto** (filtro della sezione 2.5: stima sotto 10 in valore assoluto, errore standard sotto 10). Sono dataset con 2-17 eventi e 2 occasioni: la stima di `time:trt` va a 17-37 e coincide fra le due implementazioni alle cifre stampate; l'errore standard, fra 2.900 e 11.800, differisce fino al 3,4% (dove il pacchetto ne dà uno: in 2 casi è `NA`). Con una verosimiglianza così piatta l'Hessiano numerico non è determinato.
- **GLMM-AR1, 10 soggetti, 8 occasioni, sotto H1, replicato 1.** Il fit con il trend non converge (codice `nlm` 3) in entrambe le implementazioni, quindi si passa al fit senza trend, che non converge neanche lui (codice 5) ma ha stima e SE di `time:trt` uguali alle 4-5 cifre (0,016410 e 0,14623 nel pacchetto; 0,016410 e 0,14622 nello studio). Lo studio lo usa perché chiede finiti solo stima ed errore standard del termine saggiato; il pacchetto chiede finiti tutti i coefficienti della media e qualcuno non lo è, quindi tiene il fit con il trend, il cui errore standard è `NA`: nessun risultato usabile. È una differenza di regola fra pacchetto e studio, presente già nella 0.3.0 (il motore non è cambiato, controllo 7).
- **GenBinaryY, 65 soggetti, 4 occasioni, sotto H1, replicato 2.** Il fit con il trend termina con codice 2 e pendenza ben sopra la soglia (15,9) in entrambe le implementazioni, ma in punti diversi: stima 0,00953 contro 0,00943, errore standard 0,125 contro 0,170. La log-verosimiglianza del pacchetto (-89,43) è sotto quella del fit senza trend (-87,46): è un falso massimo, quello che la condizione di coerenza intercetta. Con il default del pacchetto il fit riportato è quello senza trend (stima 0,0042, errore standard 0,313).

Fra i 119 confronti in cui lo studio ha un risultato usabile, le differenze sono quindi 2 (1,7%). L'esito del cancello resta "non passa": la tolleranza era fissata prima e non si cambia.

Descrizione preregistrata (non cancello): con il default del pacchetto (`check_coherence = TRUE`) il fit riportato cambia in 4 dei 127 dataset (GenBinaryY: 10 soggetti e 8 occasioni sotto H0 e H1, replicato 1; 25 soggetti e 4 occasioni sotto H1, replicato 3; 65 soggetti e 4 occasioni sotto H1, replicato 2); in tutti il percorso passa da "trend fit kept" a "constant fit used" e l'errore standard sale da 0,07-0,21 a 0,21-0,42 (`out/c9b_summary.txt`). Il riassunto di `c9_mmm_repro.R` stampava `NA` per questo conteggio (valori mancanti non gestiti); il numero giusto è quello di `c9b`.

### Correzioni a strumenti e test dopo il primo commit (dichiarate)

1. `test-decision.R`: il dataset di prova con 10 soggetti e 2 occasioni (seme 3) veniva rifiutato prima della scelta (esito che varia in un solo soggetto), quindi le cautele erano vuote e i controlli negativi passavano senza verificare nulla. Seme 4, con un'attesa esplicita che il piano non sia rifiutato; lo stesso nell'aiuto del test delle note di regione.
2. `test-decision.R`: il controllo negativo cercava "design points with 10 clusters", che compare anche nell'avviso del modello marginalizzato a 10 soggetti ("... of the 156 simulated design points with 10 clusters"); ora cerca "at 2 design points with 10 clusters".
3. `c5_constants.R`: il modello del testo vecchio "is generalized estimating equations with" avrebbe trovato anche la nuova descrizione dell'alternativa nel manuale (una riga sola). Ora cerca la frase della vecchia raccomandazione ("occasions, is generalized estimating", o una riga che inizia con "is generalized estimating equations with"). Corretto prima di eseguire il controllo.
4. `c9b_diagnose.R`: aggiunto dopo il fallimento del controllo 3, per spiegarlo; non cambia l'esito.

### Cosa cambia per chi usa il pacchetto

- `auto_cluster_fit()` sceglie `"mmm"` (modello marginalizzato con il controllo di convergenza, t); la ricaduta è `"gee_ar1_fg"`, con le sue cautele; se falliscono entrambe, il rifiuto dà i due motivi e suggerisce `flip_test()`.
- `cluster_strategies()$alternative` è `"gee_ar1_fg"` (con la correlazione di `fewclust_control()`).
- Nuova cautela per `"gee_ar1_fg"` con meno di 15 soggetti e 2 occasioni; nuova nota di regione (10 soggetti, 3 o più occasioni: GEE più vicina al 5%); tolta la nota "slightly closer to 5%"; la cautela sui buchi nelle serie vale anche per il modello marginalizzato.
- `flip_test()` su esito costante: `NA` con motivo "constant outcome ...".
- Conteggi delle cautele dalla tabella corretta (Tabella 5 del paper v6), generati da `dev/make_constants.R`.
- Default di `check_coherence` invariato (`TRUE`).


## Decisioni del PI (2026-11-02)

Prese nella chat dello step S9 del manoscritto, sul riassunto di quello step:

1. Le 2 differenze del controllo 3 sono accettate come spiegate; la versione 0.3.1 e' chiusa, con una sezione "Agreement with the study" in `NEWS.md`.
2. La condizione di coerenza resta attiva per default (`check_coherence = TRUE`); `NEWS.md` dice che non e' nella regola dello studio e come spegnerla.
3. Il manoscritto (v7) cita fewclust come implementazione della raccomandazione, con una frase sulla condizione di coerenza. La condizione per citarlo e' stata rifatta dal ramo del manoscritto con la 0.3.1 installata dal commit c688ecc: scelta `mmm` su 2 dataset su 2 (`R/s9/s9_fewclust_condition.R`, `tables/V6_S9_fewclust_condition.csv` su `morris-v4`).

Dopo queste decisioni cambiano solo `NEWS.md` e questo file: codice, costanti e test sono quelli del commit c688ecc, quindi i controlli non si rifanno.
