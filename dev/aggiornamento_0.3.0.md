# Aggiornamento di fewclust alla versione 0.3.0

Preparato il 2026-11-02 nella sessione parallela sul pacchetto (ramo `package`, worktree `~/MMM-package`), seguendo il brief `BRIEF_fewclust_0.3.0.md`. Il ramo del manoscritto (`morris-v4`, `~/MMM`) si legge soltanto.

**Stato di questo file.** Le sezioni da "Cosa cambia" a "Test esistenti che cambiano" sono scritte e committate *prima* di toccare il codice del pacchetto. La sezione "Esiti" si aggiunge alla fine, con i numeri letti dai file in `dev/checks_0.3.0/out/`.

## Obiettivo

Portare fewclust dalla 0.2.0 (commit `8538cad`) a una 0.3.0 che dica le stesse cose del manoscritto v5 finale (commit `4b371af` di `morris-v4`), senza cambiare i calcoli che il paper ha preso dal pacchetto (le correzioni del sandwich di `gee_vcov()`).

## Cosa cambia (piano)

1. **Strategia raccomandata.** La scelta automatica diventa la GEE con correlazione di lavoro AR(1), la correzione di Fay e Graubard e il riferimento t con N - p gradi di libertà, a ogni numero di cluster e occasioni (Figura 1, sezione 4.1 e conclusioni del manoscritto v5). Il modello marginalizzato multilivello con il controllo di convergenza e la t è l'alternativa: il piano la nomina (`alternative = "mmm"`) e la si usa se la GEE non restituisce un risultato usabile. Con 2 occasioni e almeno 15 cluster una nota dice che lì il modello marginalizzato era un po' più vicino al 5% e la GEE conservativa (sezione 3.4; decisione del PI sulla parola "comparable", D-61 nel ramo del manoscritto).
2. **Test a inversione di segno.** Nuova funzione `flip_test()`, che restituisce solo il p-value. Calcolo identico a quello del paper (sezione 2.4, supplemento S7; codice `R/s7c/s7c_tests.R`, `R/s7d/s7d_flip.R`, `R/s7d2/s7d2_ci.R` del ramo `morris-v4`): modello senza il termine saggiato (`glm.fit`, tolleranza 1e-12, al massimo 100 iterazioni), contributi efficaci allo score per soggetto con la proiezione sugli altri coefficienti, tutte le 2^N assegnazioni di segno fino a 12 soggetti, altrimenti 1.999 assegnazioni casuali più l'identità. Nessuna risposta (p-value `NA` con il motivo) se il modello ridotto non converge, diverge (un coefficiente di almeno 1.000 in valore assoluto), ha informazione singolare, o se tutti i contributi allo score sono sotto 1e-6 (regole della sezione 2.4, decisione D-59 del ramo del manoscritto). L'intervallo ottenuto invertendo il test **non** si offre (sezione 4.1). Con 2-3 occasioni e al massimo 25 cluster il piano lo indica come complemento (`complement = "flip_test"`).
3. **Correlazione di lavoro.** Il default di `fewclust_control(corstr = )` passa da `"exchangeable"` a `"ar1"`. Se l'utente sceglie un'altra correlazione con la strategia automatica, una cautela dice che la regola dello studio sceglie AR(1) e riporta il conteggio della Tabella 4 per quella correlazione.
4. **La GEE come nello studio** (trovato rileggendo il codice dello studio, non nel brief). Lo studio stima la GEE con `geepack::geeglm(..., scale.fix = TRUE)` su righe ordinate per soggetto e poi per tempo (`R/s7d/s7d_common.R::one_cor`). La 0.2.0 non fissa la scala e ordina solo per soggetto: con AR(1) l'ordine delle righe dentro il soggetto definisce i ritardi, e la scala stimata cambia la stima del parametro di correlazione. La 0.3.0 fa come lo studio. Decisione tecnica dell'agente: è l'attuazione del metodo del paper, non un metodo nuovo. Il calcolo di `gee_vcov()` non cambia.
5. **Nomi delle strategie** (necessario per la condizione pre-registrata del paper, sotto). Le strategie GEE riportano la correlazione di lavoro nel nome, come nella Tabella 4: `"gee_ar1_fg"`, `"gee_exch_md"`, `"gee_ind_kc"` e così via. In ingresso restano validi i nomi della 0.2.0 (`"gee_fg"` usa la correlazione di `fewclust_control()`).
6. **Costanti dello studio.** `R/constants.R` passa ai 2.376 punti di disegno e ai quattro meccanismi. I conteggi (punti chiaramente sopra 7,5% e sotto 2,5%, anche a 10 cluster) e le quote di risultati usabili si generano dalla Tabella 4 (`~/MMM/tables/V5_S7f_table4.csv`), per strategia e correlazione di lavoro, e non si ricopiano a mano. Le cautele con riferimento normale non citano più conteggi: la Tabella 4 del v5 non ne ha; dicono che i conteggi si riferiscono alla t. Il pacchetto non cita i risultati con N - 2 gradi di libertà, quindi `V5_S7f_rule_t2.csv` non serve.
7. **Controllo di convergenza del modello marginalizzato.** La regola resta quella del paper (soglia 1e-3, codici `nlm` diversi da 1 e 2, ricaduta alla deviazione standard latente costante, t con G - p). Testo e documentazione si allineano al paper: due modi di fallimento, la non convergenza e la convergenza con errore standard collassato; nessun "spostamento della stima" (supplemento S16, sezione 2.4). Il fit registra ciò che la specifica chiede: percorso (`"trend fit kept"`, `"constant fit used"`, `"constant fit failed, trend fit kept"`), motivi, e il fit con il trend (stima, errore standard, codice, pendenza minima) anche quando è sostituito, così l'utente può riportare quante volte il controllo è intervenuto (sezione 4.2).
8. **Errore standard del modello marginalizzato.** Il pacchetto riporta solo l'errore standard basato sul modello, come il paper; nessun errore standard robusto. La documentazione dice che nelle applicazioni cliniche del paper quell'errore standard era più piccolo della dispersione delle stime (sezione 4.4).
9. **Buchi nelle serie.** Con AR(1) i ritardi sono le osservazioni consecutive del soggetto. Una cautela segnala i soggetti con un'occasione mancante in mezzo alla serie (nello studio le serie erano complete o troncate alla fine).
10. **Documentazione.** `DESCRIPTION`, `NEWS.md` (voce 0.3.0), `README.md`, le due vignette e `inst/CITATION` riscritti sulla nuova raccomandazione, con il riferimento al manoscritto v5. Si tolgono dal ramo i 71 file `._*` di macOS sotto `fewclust/` (residui del portatile, che `R CMD check` segnala come file nascosti).

### Una scelta per il PI (aperta)

La 0.2.0 aggiunge alla regola del paper una terza condizione, attiva per default (`check_coherence = TRUE`): un fit con il trend la cui log-verosimiglianza è sotto quella del fit senza trend è trattato come non convergente. Il paper non la contiene e non l'ha valutata. La documentazione della 0.2.0 afferma che "catches the convergence with a collapsed standard error the study describes, which the two published conditions do not see": il paper dice invece che le due condizioni pubblicate coprono entrambi i modi di fallimento. Il testo si corregge comunque. Se tenerla attiva per default lo decide il PI, dopo la previsione P1 qui sotto. Finché non decide, il default della 0.2.0 resta.

## Controlli (cancelli)

Ognuno è binario; se uno fallisce, la 0.3.0 non si chiude. Script in `dev/checks_0.3.0/`, lanciati da `run_checks.sh`; risultati in `dev/checks_0.3.0/out/`. Libreria di sviluppo `~/R/fewclust_dev_lib`; `~/R/s7d_lib` non si tocca. Al massimo 4 core.

1. **Test a inversione di segno riprodotto** (`c1_study_repro.R`). Dataset dello studio rigenerati dai loro semi con il codice di `morris-v4`: scenario M1, quattro meccanismi, N in {10, 25}, T in {2, 10}, sotto H0 e H1, replicati 1-3 (96 dataset); scenari M2 e M5 con covariata, 6 dataset; e fino a 3 dataset per cella fra quelli in cui il test del paper non era calcolabile (M1, N = 10, T = 2, otto celle). Su ciascuno `flip_test(seed = 20261022 + i)` deve dare il p-value salvato (`p0_guard` in `results_v4/S7d2/rows/`) entro 1e-12, oppure nessuna risposta con lo stesso motivo. Passa se coincidono tutti, con almeno 20 dataset, e se i dati rigenerati hanno lo stesso numero di righe e di eventi del run salvato.
2. **Strategia GEE riprodotta** (stesso script). Sugli stessi dataset, dove il pacchetto non rifiuta e lo studio ha una stima GEE AR(1) con errore standard di Fay e Graubard finito, `auto_cluster_fit()` deve scegliere `"gee_ar1_fg"` e dare stima ed errore standard del termine `time:trt` uguali a quelli salvati (`ar1_est`, `ar1_se_fg`), con tolleranza relativa 1e-8. Passa se coincidono tutti, su almeno 20 dataset confrontati. I rifiuti del pacchetto e i dataset dove la GEE dello studio non era usabile si contano e si descrivono.
3. **Regola come nel paper** (`c2_rule_fig1.R` e `tests/testthat/test-decision.R`). Per N in {10, 25, 145} e T in {2, 3, 10}, `cluster_strategies()` deve dare: strategia `"gee_ar1_fg"`, riferimento t, alternativa `"mmm"`, complemento `"flip_test"` se e solo se T ≤ 3 e N ≤ 25, nota "slightly closer to 5%" per il modello marginalizzato se e solo se T = 2 e N ≥ 15. Passa se tutte e 9 le celle coincidono con la Figura 1.
4. **`R CMD check --as-cran`** (`c7_rcmdcheck.sh`): nessun errore e nessun warning; ogni NOTE motivato in `cran-comments.md`.
5. **Test del pacchetto** (`testthat`): tutti passano. I test delle correzioni (`test-corrections.R`) e del motore del modello marginalizzato (`test-mmm.R`, salvo la richiesta esplicita `strategy = "mmm"` dove prima si passava dalla scelta automatica) restano con le stesse attese numeriche. Il nuovo `test-mmm-check.R` verifica sulle sei fixture del controllo di convergenza (`dev/vedovelli_fixtures/`) ramo, stima ed errore standard usati e p-value con la t su 21 gradi di libertà, con la regola del paper (`check_coherence = FALSE`); tolleranze: stima 1e-6, errore standard 1e-5, p-value 1e-5.
6. **I calcoli presi dal paper non cambiano** (`c3_invariance_*.R`). `gee_vcov()` (non corretta, Fay e Graubard, Mancl e DeRouen, Kauermann e Carroll) su 100 dataset simulati, con le tre correlazioni di lavoro e con `scale.fix` vero e falso, e il motore del modello marginalizzato con un controllo esplicito su 40 dataset, devono dare risultati identici bit per bit fra la 0.2.0 (`~/R/s7d_lib`) e la 0.3.0.
7. **Condizione pre-registrata del paper per citare il pacchetto** (`c4_paper_condition.R`; decisioni D-43 e D-65 del ramo del manoscritto). Lo stesso codice di `R/s7f/s7f_fewclust_condition.R`, con la 0.3.0, sugli stessi due dataset: `condition_met` vero su entrambi. Lo script scrive solo nella cartella dei controlli del pacchetto, non in `~/MMM/tables`.
8. **Costanti dalla Tabella 4** (`c5_constants.R`). Ogni conteggio e quota che il pacchetto cita è uguale al valore di `V5_S7f_table4.csv`; nessun testo del pacchetto cita ancora i 1.896 punti del v4; `DESCRIPTION` ha la versione 0.3.0 e `NEWS.md` si apre con la 0.3.0.

Integrità: `run_checks.sh` verifica che commit e stato di `~/MMM` non cambino durante i controlli; `c1_study_repro.R` confronta l'md5 dei file dei risultati letti prima e dopo.

## Previsione P1: quanto pesa la condizione di coerenza

- **Quantità**: la quota di dataset in cui la condizione di coerenza cambia il fit riportato dal modello marginalizzato, rispetto alla regola del paper.
- **Dataset**: scenario M3 (deviazione standard latente che cresce nel tempo), quattro meccanismi, N in {25, 65}, T in {4, 8}, sotto H0, replicati 1-50: 800 dataset rigenerati dai semi (`c6_coherence.R`).
- **Valore atteso**: circa 3%, come dice la documentazione della 0.2.0.
- **Intervallo di rilevanza**: da 1% a 6%. Fuori da questo intervallo la frase "about 3 percent" va cambiata. Con 800 dataset e una quota del 3%, 2 MCSE valgono circa 1,2 punti, quindi l'esito non è indeterminato per costruzione.
- Si riportano anche, fra i fit cambiati, la differenza della stima di `time:trt` in unità di errore standard e il rapporto fra gli errori standard: servono al PI per decidere il default.

## Test esistenti che cambiano

Le attese che codificano la vecchia raccomandazione cambiano nel primo commit, insieme ai controlli:

- `test-decision.R`: scelta automatica `"gee_ar1_fg"` con alternativa e complemento; la griglia 3 × 3 della Figura 1; conteggi su 2.376 punti, per correlazione di lavoro; cautele con riferimento normale senza conteggi; la cautela del modello marginalizzato a 10 cluster con 18 punti su 156 (Tabella 4); i buchi nelle serie con AR(1).
- `test-fit.R`: la scelta automatica è una GEE AR(1) con Fay e Graubard; la ricaduta va dalla GEE al modello marginalizzato (prima il contrario); se falliscono entrambe, il rifiuto suggerisce `flip_test()`.
- `test-compare.R`, `test-methods.R`: nomi delle strategie con la correlazione di lavoro; la scelta del modello marginalizzato ora è un override.
- `test-control.R`: default `corstr = "ar1"`.
- `test-mmm.R`: dove si passava dalla scelta automatica si chiede `strategy = "mmm"`; nessuna attesa numerica cambia.
- Nuovi: `test-flip.R` (confronto con una copia del codice dello studio, enumerazione, segni casuali, regole di calcolabilità, ordine delle righe, seme), `test-mmm-check.R` (fixture del controllo di convergenza), in `test-gee.R` la GEE come nello studio.

## Esiti

Scritti il 2026-11-02 dopo l'ultimo giro completo dei controlli sul codice finale (`run_checks.sh all`, log in `out/run_checks.log`). Ogni numero è letto dai file di `dev/checks_0.3.0/out/`.

### In breve: il pacchetto ora dice le stesse cose del paper?

**A cosa serviva.** Il pacchetto fewclust raccomandava ancora il modello marginalizzato multilivello come prima scelta, mentre il paper finale raccomanda la GEE con correlazione di lavoro AR(1), la correzione di Fay e Graubard e la t. Per questo il paper non poteva citarlo come implementazione della sua raccomandazione. L'aggiornamento allinea il pacchetto al paper senza cambiare i calcoli che il paper ha preso dal pacchetto.

**Cosa abbiamo trovato.**

1. **Il pacchetto ora sceglie la strategia del paper, e la applica come nello studio.**
   - A ogni numero di soggetti e di visite sceglie la GEE con correlazione AR(1), Fay e Graubard e la t. Indica il modello marginalizzato come alternativa e lo usa se la GEE non dà un risultato. Con 2-3 visite e fino a 25 soggetti propone in più il test a inversione di segno.
   - Su 9 combinazioni di soggetti e visite la scelta è quella della Figura 1 in tutte e 9.
   - Su 51 dataset dello studio rigenerati dai semi, stima ed errore standard sono identici a quelli salvati dallo studio, fino all'ultima cifra.
   - Il controllo che il paper aveva fissato per poter citare il pacchetto, fallito con la versione precedente, ora passa su entrambi i dataset di prova.
2. **Il test a inversione di segno è nel pacchetto e dà gli stessi p-value del paper.** Su 126 dataset dello studio il p-value coincide esattamente con quello salvato, oppure entrambi dicono che il test non è calcolabile, per lo stesso motivo. Il pacchetto restituisce solo il p-value e non propone l'intervallo, come il paper.
3. **I numeri del paper non cambiano.** Le correzioni del sandwich calcolate dal pacchetto sono identiche, bit per bit, a quelle della versione usata per il paper, su 600 combinazioni di dataset e correlazioni di lavoro; lo stesso vale per il modello marginalizzato su 40 dataset. Il controllo formale per la sottomissione al repository dei pacchetti R non dà né errori né avvisi.
4. **Una differenza voluta per chi usa il pacchetto: la GEE ora è stimata come nello studio.** Lo studio fissava il parametro di scala e ordinava le visite di ogni soggetto nel tempo, mentre la versione precedente del pacchetto no. Con la correlazione AR(1) questo cambia un poco stime ed errori standard rispetto alla versione precedente.
5. **Il controllo di convergenza del modello marginalizzato ha, nel pacchetto, una condizione in più che il paper non ha.** È la condizione di coerenza: un fit con il trend la cui verosimiglianza è sotto quella del fit senza trend si scarta.
   - Su 796 dataset dello studio con la deviazione standard latente che cresce nel tempo, cambia il risultato in 13, cioè circa 2 casi su 100.
   - In 12 di questi casi il fit scartato aveva un errore standard fra un terzo e tre quinti di quello del fit senza trend. È proprio il difetto che il controllo deve intercettare, e che la regola del paper in questi casi lascia passare.

**Cosa devi decidere.**

1. **Teniamo attiva per default la condizione in più del controllo di convergenza? Consiglio di tenerla.**
   - Ragioni: nei pochi casi in cui interviene corregge un errore standard troppo piccolo. Inoltre non riguarda la strategia raccomandata, che è la GEE.
   - Se la teniamo: il pacchetto lo dichiara come aggiunta propria, e chi vuole la regola esatta del paper la spegne con un'opzione. Il modello marginalizzato del pacchetto può allora differire da quello del paper in circa 2 dataset su 100.
   - Se la spegniamo: il modello marginalizzato del pacchetto è esattamente quello del paper, compresi quei casi con errore standard collassato.

**Comunque.** Quando nella chat del manoscritto installeremo la nuova versione e rifaremo la pipeline, la tabella che registra la strategia scelta dal pacchetto cambierà, come atteso: ora sceglie la GEE AR(1) con Fay e Graubard. Le altre cifre del paper non devono cambiare. Passo a quella chat anche un'osservazione che riguarda il paper e non il pacchetto. Fra i dataset dello studio, 9.646 hanno l'esito tutto 0 o tutto 1 e contano come risultato usabile della GEE: la stima che ne esce è solo errore di arrotondamento. Togliendoli, la quota di risultati usabili della GEE passa dal 92,1% al 91,9% su tutto lo studio e dal 45,9% al 43,7% con 2-3 visite e fino a 25 soggetti.

Appena decidi questo punto, la 0.3.0 è chiusa.

### Cancelli

| Controllo | Criterio | Osservato | Esito |
|---|---|---|---|
| 1. Test a inversione di segno riprodotto | tutti i p-value entro 1e-12 o non calcolabili con lo stesso motivo; ≥ 20 dataset; dati identici al run salvato | 126 su 126 (73 con tutte le 2^N assegnazioni, 51 con 1.999 assegnazioni casuali, 49 non calcolabili con lo stesso motivo: 43 "score numerically zero", 6 "zero score variance"; 2 non generabili, qui come nel run salvato); differenza massima 0; md5 dei file letti invariato | passa |
| 2. Strategia GEE riprodotta | stima e SE di `time:trt` uguali con tolleranza relativa 1e-8; ≥ 20 confronti | 51 su 51, differenza massima 0 | passa |
| 3. Regola come nella Figura 1 | 9 celle su 9 | 9 su 9 | passa |
| 4. `R CMD check --as-cran` | 0 errori, 0 warning | `Status: 1 NOTE` (validazione HTML del manuale saltata: `tidy` non installato) | passa |
| 5. Test del pacchetto | tutti passano | 15 file, 632 attese, 0 fallite, 0 errori | passa |
| 6. Calcoli presi dal paper invariati | identici bit per bit fra 0.2.0 e 0.3.0 | `gee_vcov()`: 600 insiemi su 600 (64 contengono un messaggio d'errore, identico anch'esso); motore del modello marginalizzato: 40 su 40 | passa |
| 7. Condizione pre-registrata del paper | `condition_met` vero sui due dataset | vero su 2 su 2 (`chosen = gee_ar1_fg`) | passa |
| 8. Costanti dalla Tabella 4 | ogni valore uguale; nessun riferimento ai 1.896 punti | 133 valori su 133 per 19 strategie; 0 righe obsolete; versione 0.3.0 | passa |

Fonti: `out/c1_summary.txt` (controlli 1 e 2), `out/c2_summary.txt`, `out/c7_summary.txt` e `out/c7_check.log`, `out/tests.log`, `out/c3_summary.txt`, `out/c4_summary.txt`, `out/c5_summary.txt`. Integrità: commit e stato di `~/MMM` invariati durante i controlli (`out/run_checks.log`).

Dettaglio del controllo 2. Dei 126 dataset, 75 non sono confrontati: 63 perché il pacchetto rifiuta prima di stimare o non ottiene un fit (24 esito che varia in un solo soggetto, 7 esito costante, 16 separazione, 6 disegno degenere a rango non pieno, 10 GEE e modello marginalizzato entrambi senza risultato), 10 perché nello studio la GEE AR(1) con Fay e Graubard non aveva un errore standard finito (lì il pacchetto passa al modello marginalizzato, come previsto), 2 perché il dataset non si genera. I rifiuti prima della stima sono regole del pacchetto (dalla 0.1.0), non del paper.

### Previsione P1: indeterminata

13 fit cambiati su 796 con entrambi i fit disponibili (4 dataset su 800 senza uno dei due fit): 1,63%, MCSE 0,45, stima ± 2 MCSE da 0,73% a 2,53%. L'intervallo è a cavallo dell'estremo inferiore (1%) dell'intervallo di rilevanza, quindi l'esito è indeterminato. Tocca solo la frase della documentazione, che ora riporta il numero misurato (13 su 796) invece di "about 3 percent". Nei 13 casi il fit con il trend era convergente, con pendenza sopra soglia. Il rapporto fra l'errore standard senza trend e quello con il trend ha mediana 2,20 (da 0,99 a 2,98; da 1,68 in su in 12 casi su 13). La differenza della stima, in unità di errore standard, ha mediana 0,005 (da -2,61 a 1,46). Fonte: `out/c6_summary.txt`, `out/c6_coherence.csv`.

**Correzione al testo pre-registrato.** Nella sezione "Una scelta per il PI" avevo scritto che il paper afferma che le due condizioni coprono entrambi i modi di fallimento. Il paper dice che i due modi esistono e che il controllo elimina l'eccesso di rifiuti (sezioni 2.4 e 3.2, supplemento S16), non che ogni errore standard collassato venga intercettato. P1 mostra 13 casi su 796 in cui non lo è. La frase della 0.2.0 che la condizione di coerenza "catches the convergence with a collapsed standard error ... which the two published conditions do not see" è dunque in parte confermata. Nella 0.3.0 è sostituita dal numero misurato.

### Correzioni a strumenti e test dopo il primo commit (dichiarate)

1. `test-mmm-check.R`: il confronto del fit con il trend salvato dallo studio si fa solo dove `nlm` ha converso (codici 1-2). Nella fixture `case_04` (codice 3, non convergente) il punto d'arresto dipende dall'implementazione: stima 0,0054328 qui contro 0,0054321 in MMLB, errore standard 0,086 contro 0,099. Ramo, stima ed errore standard usati e p-value, cioè ciò che il cancello 5 chiedeva, coincidono in tutte e sei le fixture.
2. `test-mmm.R`: il confronto fra tempo in anni e in giorni si fa su ramo, esito e percorso del controllo, non sull'intero elenco. L'elenco ora contiene anche il fit con il trend, che con il tempo in giorni segue un altro percorso dell'ottimizzatore.
3. `c1_study_repro.R` e `c6_coherence.R`: sotto R 4.6 il `library(fewclust)` del codice dello studio sostituiva la 0.3.0 caricata con la 0.2.0 di `~/R/s7d_lib`. Gli script ora la ricaricano dalla libreria di sviluppo e verificano versione e percorso.
4. `c1_study_repro.R`: un dataset che il codice dello studio non riesce a generare (2 casi, GLMM-AR1, 10 soggetti, 2 visite, replicato 2) si confronta con il flag salvato (`generated = FALSE` anche nel run salvato). Un errore di `auto_cluster_fit()` si registra come mancato confronto.
5. `c7_rcmdcheck.sh`: `~/bin/tidy` è un wrapper rotto, che punta a un file temporaneo di un'altra sessione e interrompeva il controllo del manuale HTML. Il comando `tidy` di R è puntato altrove (`R_TIDYCMD`), per cui la validazione HTML è saltata: è l'unico NOTE.
6. Pacchetto, trovato dal cancello 1: `flip_test()` si fermava con un errore sui disegni degeneri, dove nessun soggetto trattato è osservato dopo la prima occasione e la colonna dell'interazione è tutta zero. Ora risponde "non calcolabile" con il motivo dello studio ("zero score variance"), come il paper. Sono 6 dei 126 dataset, esclusi da ogni denominatore del paper. Test aggiunto in `test-flip.R`.

### Osservazione per il manoscritto (non per il pacchetto)

Calcolo in sola lettura sui replicati salvati (`c8_constant_outcome.R`, `out/c8_summary.txt`). I replicati con disegno stimabile sono 4.749.358. Fra questi, 9.711 hanno l'esito costante e 9.646 di questi contano come usabili per la GEE AR(1) con Fay e Graubard secondo il filtro della sezione 2.5: la stima è dell'ordine di 1e-14, l'errore standard è finito. La quota media di risultati usabili ricalcolata coincide con la Tabella 4: 92,121% su 2.376 punti e 45,944% con 2-3 visite e al massimo 25 cluster. Senza quei replicati scende a 91,897% e a 43,72%. La quota più alta di replicati con esito costante in un punto di disegno è 13,13% (Copula-AR1, M1, 10 cluster, 2 visite). Il test a inversione di segno li esclude già per regola (decisione D-59 del ramo del manoscritto); per la GEE e, forse, per il modello marginalizzato no. Va valutato nella chat del manoscritto: lì si decide se cambia qualcosa (con 2-3 visite e fino a 25 soggetti la frase "46%" diventerebbe "44%"; la scelta della regola non cambia, perché la differenza fra le due strategie resta entro il margine).

### Cosa cambia per chi usa il pacchetto

- `auto_cluster_fit()` sceglie `"gee_ar1_fg"`; la ricaduta è `"mmm"`; se falliscono entrambe, il rifiuto suggerisce `flip_test()`.
- `fewclust_control(corstr = "ar1")` è il nuovo default; le strategie GEE si chiamano con la correlazione di lavoro (`"gee_ar1_fg"`, `"gee_exch_md"`, ...), e i nomi vecchi restano validi in ingresso.
- La GEE è stimata con `scale.fix = TRUE` e righe ordinate per soggetto e tempo: stime ed errori standard delle strategie GEE possono differire da quelli della 0.2.0.
- Nuova `flip_test()` (solo p-value). Il piano ha `$alternative` e `$complement`. Il fit del modello marginalizzato registra percorso, motivi e i due fit del controllo di convergenza.
- Cautele con i conteggi della Tabella 4 del paper finale, per strategia e correlazione di lavoro.

### Verifica finale (nella chat del manoscritto)

Installare la 0.3.0 da questo commit in `~/R/s7d_lib`, rieseguire la pipeline del paper e confrontare il registro dei numeri con l'attuale: 0 differenze attese. Cambieranno invece `tables/V5_S7f_fewclust_condition.csv` (`chosen = gee_ar1_fg`, `condition_met = TRUE`), che non è nel registro, e la riga di `R/s7f/s7f_inbreve_sources.R` che si aspetta `"mmm"`. Solo allora si aggiornano versione e commit nella Tabella S2 e in "Code and data availability".
