# fewclust 0.4.1: il fit senza pendenza conta solo se converge (D-90)

Preparato il 2026-11-03. Il PI ha deciso (D-90, `~/MMM/docs/decisions/D-90.md`, che sostituisce D-32 nel punto 3) che nel controllo di convergenza del modello marginalizzato multilivello un fit "al bordo" si sostituisce con il fit senza la pendenza della deviazione standard latente **solo se** `nlm` dichiara la convergenza di questo secondo fit (codice 1 o 2) **e** stima ed errore standard sono finiti; altrimenti il modello marginalizzato non dà un risultato usabile per quel dataset. Si parte dalla 0.4.0 (commit `e7bb055`, ramo `package`). Il ramo del manoscritto (`~/MMM`, `/mnt/micron/MMM`, ramo `morris-v4`) e la libreria congelata dello studio (`~/R/s7d_lib`) non si toccano.

**Stato di questo file.** Le sezioni fino a "Test" sono scritte e committate *prima* di qualsiasi modifica del codice. Gli esiti si aggiungono alla fine, senza cambiare cancelli e attese.

## Cosa cambia

1. **`fit_mmm()` (R/mmm.R), i due rami del controllo.** Oggi un fit al bordo (codice `nlm` diverso da 1 e 2, o pendenza della deviazione standard latente sotto `check_threshold`, o, con `check_coherence = TRUE`, log-verosimiglianza sotto quella del fit senza pendenza) è sostituito dal fit senza pendenza quando questo ha stime ed errori standard finiti, *con qualunque codice `nlm`*; se il fit senza pendenza non è finito si tiene il fit al bordo. Dalla 0.4.1:
   - fit senza pendenza con codice 1 o 2 e stime ed errori standard finiti: si usa (percorso `"constant fit used"`, come prima);
   - fit senza pendenza con codice 3, 4 o 5: **nessun risultato** (percorso nuovo `"no result"`, `refit_status = "not converged"`, il motivo dice il codice);
   - fit senza pendenza non finito o fallito con errore: **nessun risultato** (`"no result"`, `refit_status = "not finite"`);
   - il fit al bordo non si riporta più in nessun caso; il percorso `"constant fit failed, trend fit kept"` e l'etichetta `"time trend (boundary fit kept)"` spariscono.
   La definizione di "al bordo" non cambia (anche la condizione di coerenza resta com'è, calcolata come nella 0.4.0 sul fit senza pendenza finito).
2. **Nessun risultato = fallimento del modello marginalizzato**, della stessa classe interna (`fc_failure`) dei fallimenti già esistenti, con il record del controllo allegato. Quindi: con `method = "auto"` si ricade sulla GEE (AR(1) e Fay-Graubard di default) con il messaggio già esistente; con `method = "mmm"` errore `fewclust_error` come per ogni altro fallimento del modello marginalizzato. Il testo (cli) dice perché: "the refit without the time trend did not converge (nlm code k)" oppure "returned no finite estimate".
3. **Stampa e sommario** (R/methods.R): l'etichetta del controllo dice "intervened (trend fit at the boundary)" solo per un fit senza pendenza convergente; per il percorso nuovo dice che il modello marginalizzato non ha risultato e perché. Con la ricaduta sulla GEE, `fit$convergence` conserva il record del controllo del modello marginalizzato e `summary()` lo mostra.
4. **Documentazione**: `fewclust_control(check = )`, `fit_fewclust()` (sezione Methods), commento di `fit_mmm()`, vignetta "Getting started", NEWS (sezione 0.4.1, con il comportamento precedente), DESCRIPTION `Version: 0.4.1`. La nota "the reported fit's optimiser returned nlm code" diventa impossibile e si toglie (anche dal filtro di `fc_fit_notes()`).
5. **Invariati**: motore (`fc_mmm_fit_once()`, verosimiglianza, Hessiana, valori iniziali), GEE e correzioni, `gee_vcov()`, `flip_test()`, `simulate_clustered_binary()`, `R/constants.R` (le costanti restano quelle delle tabelle dello studio di S9; si aggiorneranno quando lo studio avrà ricalcolato le tabelle con D-90, non in questo passo).

## Controlli (cancelli)

Script in `dev/checks_0.4.1/`, lanciati da `run_checks.sh`; uscite piccole in `dev/checks_0.4.1/out/`, oggetti intermedi in `/home/user/.claude-science-scratch/fc041/` (fuori dal repo). Riferimento: la 0.4.0 installata da `e7bb055` (git archive) in una libreria temporanea dello scratch; la 0.4.1 in una libreria nuova `~/R/fewclust_041_lib`; dipendenze lette dalla libreria renv del progetto (sola lettura). Al massimo 6 processi R contemporanei (host condiviso, limite termico).

**A. Stessi numeri della 0.4.0 fuori dai dataset previsti** (`c1_identity_*.R`, `c2_study_*.R`). Sui 60 dataset del controllo 1 della 0.4.0 (rigenerati con `simulate_clustered_binary()` e gli stessi semi; 14 stime per dataset: `auto`, `mmm`, GEE 3 correlazioni × 4 correzioni; controllo di default) e sui dataset dello studio del controllo C (sotto, `mmm` e `auto` con `check_coherence = FALSE`): l'insieme dei casi che cambiano è *esattamente* quello previsto dal record della 0.4.0 (percorso `"constant fit used"` con fit senza pendenza di codice diverso da 1 e 2, oppure `"constant fit failed, trend fit kept"`). Fuori da quell'insieme stima, errore standard, statistica, gradi di libertà e p-value sono `identical()` alla 0.4.0, e così il metodo usato. Dentro: `mmm` dà un errore `fewclust_error` il cui testo nomina il motivo (codice `nlm` o non finito); `auto` dà la GEE AR(1) con Fay-Graubard `identical()` a `method = "gee"` della 0.4.1 sullo stesso dataset. `flip_test()` e `gee_vcov()` identici sui 60 dataset; md5 di `R/constants.R` identico alla 0.4.0. I fit registrati (`trend`, `constant`: codice, stime, errori standard, log-verosimiglianza) sono identici alla 0.4.0 su tutti i dataset (il motore non cambia). Passa se tutto coincide.
**B. La regola applicata è D-90, sui fit del pacchetto.** Su ogni stima `mmm`/`auto` dei controlli A e C, l'esito della 0.4.1 coincide con la regola D-90 applicata ai fit registrati dal pacchetto stesso: al bordo se codice del fit con pendenza ∉ {1, 2} o pendenza < 1e-3 (o, col default, la condizione di coerenza); se al bordo, fit senza pendenza se codice ∈ {1, 2} e finito, altrimenti nessun risultato. Passa se 0 discordanze.
**C. Dati e integrità.** Campione di dataset dello studio (sotto): ogni dataset rigenerato con i generatori e i semi dello studio (`R/s7d2/s7d2_common.R` di `morris-v4`, in sola lettura, come in `dev/checks_0.3.1/c9_mmm_repro.R`) ha `n_obs` e `sum_y` uguali a quelli salvati in `results_v4/S7d2/rows`; md5 dei file dello studio letti uguale prima e dopo; nessun file scritto fuori da `dev/checks_0.4.1/`, dallo scratch e da `~/R/fewclust_041_lib`. Passa se tutto è vero.
**D. Pacchetto.** Test: tutti passano (0 fallite, 0 errori), compresi i test nuovi di `test-mmm-check.R` che forzano ogni ramo (fit senza pendenza convergente → usato; codice 3-5 → nessun risultato; non finito → nessun risultato; `auto` → GEE). `lintr::lint_package()` come nella 0.4.0: 0 lint in `R/` (escluso `R/constants.R`). `R CMD check --as-cran`: 0 errori, 0 warning (le NOTE si riportano). Passa se tutto è vero.

### Campione dello studio (controllo C)

Dal file delle chiavi `~/MMM/results_v4/S11/cv_keys.rds` (dataset in cui il modello marginalizzato con il controllo non dà risultato secondo D-90) e dai fit dello studio (`results_v4/S6cc/cc_all.rds` per i tre meccanismi della griglia principale, `results_v4/S7e/rows/*.rds$cc` per la copula), scenari M3-M5. Seme di campionamento 20261103, campionamento casuale semplice dentro gli strati:

- **chiavi**: strati meccanismo (4) × motivo (`refit not converged`, `refit not finite`) × scenario; 30 per strato in M3 (240), 10 per strato in M4 e M5 (80 + 80);
- **al bordo con fit senza pendenza convergente** (non chiavi: fit A al bordo, fit B finito con codice 1 o 2, `used == "B"`): strati meccanismo × scenario; 60 per strato in M3 (240), 20 in M4 e M5 (80 + 80).

Totale 800 dataset (480 in M3). Ogni dataset si stima con la 0.4.1 (e con la 0.4.0 per il cancello A) con `binY ~ time * trt` (M5: `+ dose`, la media dello studio), `id`, `time`, `method = "mmm"` e `"auto"`, `fewclust_control(check_coherence = FALSE)` (la regola dello studio).

**Limite dichiarato prima del calcolo.** Il pacchetto stima un solo modello di dipendenza: deviazione standard latente con pendenza sull'indice dell'occasione, transizione costante. È il modello dello studio in M3. In M4 (transizione `~ trt`) e M5 (deviazione standard latente `~ time + dose`, transizione `~ trt + steroids`) lo studio stima un modello che il pacchetto non può stimare: il confronto con lo studio in M4-M5 è descrittivo (quanti "nessun risultato", quante concordanze), non una verifica. Le attese di concordanza valgono per M3.

## Attese (previsioni)

Tolleranze del controllo 3 della 0.3.1: stima `all.equal(tolerance = 1e-6)`, errore standard `all.equal(tolerance = 1e-5)`.

1. **Chiavi M3 (240)**: la 0.4.1 con `method = "mmm"` non dà risultato (errore `fewclust_error`) in tutti, tranne i dataset in cui l'ottimizzatore del pacchetto arriva a un codice diverso da quello di `MMLB::mm` (fit con pendenza non al bordo per il pacchetto, oppure fit senza pendenza convergente e finito per il pacchetto). Attesa: almeno 95% senza risultato; ogni eccezione elencata con i codici dello studio e del pacchetto. Con `method = "auto"` ogni dataset senza risultato dà la GEE.
2. **Al bordo con fit senza pendenza convergente, M3 (240)**: la 0.4.1 usa il fit senza pendenza e stima ed errore standard di `time:trt` concordano con il fit B dello studio entro le tolleranze in almeno 95%; ogni eccezione elencata con i codici.
3. **M4-M5 (320)**: solo descrittivo (vedi il limite).
4. **Confronto di S9 ripetuto con la 0.4.1** (`c3_coherence.R`: stesso insieme di 162 dataset M3 di `dev/checks_0.3.1/c9_mmm_repro.R`, stesse definizioni di `c9b_diagnose.R`, chiamata `fit_fewclust(..., method = "auto")`): dataset confrontati 127 (come nella 0.3.1, salvo chiavi in cui anche la GEE fallisce) e dataset in cui la condizione di coerenza cambia il fit riportato 4 (come nella 0.3.1: dove la coerenza interviene con un fit senza pendenza non convergente la 0.4.1 ricade sulla GEE, e il fit riportato resta diverso da quello senza la condizione). File `dev/checks_0.4.1/out/c3_coherence.csv` con le colonne di `tables/V6_S9_fewclust_coherence.csv` (`item`, `value`, `source`), stesse voci, più le stesse voci calcolate contro il fit dello studio secondo D-90 (suffisso `_d90`).

## Test

In `tests/testthat/test-mmm-check.R`: si toglie il test "a failing constant fit keeps the trend fit and says so" (comportamento abolito) e si aggiungono, sui dataset `fixtures/vedovelli` e con `testthat::local_mocked_bindings()` su `fc_mmm_fit_once()` per forzare il fit senza pendenza:

1. fit senza pendenza convergente → usato (percorso `"constant fit used"`, `refit_status = "converged"`, stessi numeri del test esistente);
2. fit senza pendenza con codice 3, 4 e 5 → `fc_failure` con percorso `"no result"`, `refit_status = "not converged"`, motivo con il codice; con `method = "mmm"` errore `fewclust_error` che nomina il codice;
3. fit senza pendenza non finito (stima `NA`) e fallito con errore → `"no result"`, `refit_status = "not finite"`;
4. `method = "auto"` → GEE identica a `method = "gee"`, messaggio `fewclust_message` con il motivo, `fit$convergence$path == "no result"`, `summary()` stampa il motivo.

Gli altri test non cambiano, salvo le attese sui testi abrogati.

## Esiti

Eseguiti il 2026-11-03 con `bash dev/checks_0.4.1/run_checks.sh` sul server (R 4.6.1, al massimo 6 processi R); uscite in `dev/checks_0.4.1/out/`. Il campione dello studio è in `out/c2_sample.csv` (800 dataset, estratto con `c2_study_sample.R`, seme 20261103), le righe confrontate in `out/c2_study_rows.csv`, le eccezioni in `out/c2_exceptions.csv`.

| Cancello | Esito | Dettaglio |
|---|---|---|
| A. Stessi numeri della 0.4.0 fuori dai dataset previsti | PASS | 60 dataset del pacchetto: 960 stime, 24 previste in cambiamento (6 dataset: 1, 25, 41, 50, 54, 57; tutti con fit senza pendenza di codice 5; `auto` e `mmm`, con e senza condizione di coerenza), 936 su 936 identiche alla 0.4.0, 24 su 24 cambiate come atteso; `flip_test()` 60/60 e `gee_vcov()` 60/60 identici; md5 di `R/constants.R` identico; fit registrati identici 228/228. Campione dello studio: 800 su 800 (221 dataset previsti in cambiamento con `mmm`, 267 a livello del motore), fit registrati identici 752/752 |
| B. La regola applicata è D-90 sui fit del pacchetto | PASS | 960/960 stime sui 60 dataset, 800/800 dataset dello studio (`mmm`, `auto`, motore) |
| C. Dati e integrità | PASS | 800/800 dataset rigenerati con `n_obs` e `sum_y` uguali a `results_v4/S7d2/rows`; md5 dei file dello studio letti uguale prima e dopo; scritture solo in `dev/checks_0.4.1/`, nello scratch e in `~/R/fewclust_041_lib` |
| D. Pacchetto | PASS | Test: 13 file, 569 attese, 0 fallite, 0 errori, 0 avvisi. Lint: 0 in `R/`. `R CMD check --as-cran`: 0 errori, 0 warning, 1 NOTE (validazione HTML del manuale saltata: HTML Tidy non installato sul server, come nella 0.4.0) |

### Previsioni

1. **Chiavi M3: confermata.** 238 su 240 senza risultato con `method = "mmm"`: 191 per il controllo di convergenza (percorso `"no result"`), 47 perché `fit_fewclust()` rifiuta il disegno prima di stimare (quasi separazione, esito che varia in un solo gruppo o in nessuno; lo studio invece li stima, D-89); sul motore da solo (`fc_parse()` + `fit_mmm()`, senza il controllo del disegno) 238 su 240 senza risultato. Le 2 eccezioni sono dati in cui l'ottimizzatore del pacchetto arriva a convergenza dove `MMLB::mm` no, quindi il controllo non interviene:
   - GenBinaryY, N = 25, T = 3, H0 (`manual`), replicato 425: studio fit A codice 3, fit B codice 3 con errore standard non finito; pacchetto fit con pendenza codice 2, pendenza 5,6 (sopra la soglia), fit tenuto;
   - GenBinaryY, N = 65, T = 2, H0, replicato 602: studio fit A codice 3, fit B codice 3 non finito; pacchetto fit con pendenza codice 2, pendenza 3,1, fit tenuto.
   Con `method = "auto"` ogni dataset senza risultato del modello marginalizzato dà la GEE, oppure, dove anche la GEE fallisce, l'errore "Neither model gave a usable estimate" (cancello A).
2. **Fit senza pendenza convergente, M3: non raggiunta come scritta.** Concordano con il fit B dello studio (stima 1e-6, errore standard 1e-5) 224 su 240 (93,3%, attesa almeno 95%). Le 16 discordanze, una per una in `out/c2_exceptions.csv`:
   - 7 rifiuti del disegno di `fit_fewclust()` (quasi separazione o esito che varia in un solo gruppo), che il piano non aveva previsto; sul motore 2 di questi concordano con lo studio, 2 non hanno risultato (fit senza pendenza con errore standard non finito nel pacchetto), 3 sono rifiutati già da `fc_parse()`;
   - 3 senza risultato per il controllo: il fit senza pendenza del pacchetto (codice 1) ha errore standard non finito, quello di `MMLB::mm` no (Hessiana numerica diversa);
   - 3 fit con pendenza tenuto: il pacchetto arriva a codice 2 dove `MMLB::mm` si ferma a codice 3, quindi per il pacchetto il fit non è al bordo;
   - 3 differenze numeriche dell'errore standard con la stessa stima: Copula-AR1 N = 20, T = 3, replicato 956 (0,109 contro 0,206; fit con pendenza codice 5 in entrambi) e due stime con errore standard di migliaia (7.728 contro 7.744; 6.682 contro 5.886), comunque non usabili.
   Senza i 7 rifiuti del disegno la concordanza è 224 su 233 (96,1%). Nessuna discordanza viene dalla regola: in ognuna il pacchetto applica D-90 ai propri fit (cancello B).
3. **M4-M5 (descrittivo).** Chiavi: senza risultato 36 su 80 in M4 e 40 su 80 in M5 (22 e 21 per il controllo); fit senza pendenza convergente: 0 su 160 concordano con il fit B. Come dichiarato prima del calcolo, il pacchetto stima lì un modello di dipendenza diverso da quello dello studio, quindi questi numeri non verificano nulla.
4. **Confronto di S9 con la 0.4.1: confermata.** `out/c3_coherence.csv`: 125 dataset confrontati (127 nella 0.3.1) e 4 in cui la condizione di coerenza cambia il fit riportato (4 nella 0.3.1, gli stessi 4 dataset). I 2 dataset in meno sono GenBinaryY, N = 65, T = 2, replicato 1, H0 e H1: lo studio sotto D-90 non ha risultato, il modello marginalizzato del pacchetto neanche e la GEE fallisce, quindi `fit_fewclust()` dà un errore e il dataset esce dal confronto (era il caso previsto dal piano). Contro il fit dello studio sotto D-90: 114 confrontati, 7 discordanze (le 7 della 0.3.1 con fit A; la discordanza GLMM-AR1 N = 10, T = 8 della 0.3.1 sparisce: ora né lo studio né il pacchetto danno un risultato); gli 11 dataset senza risultato per lo studio sono senza modello marginalizzato anche nel pacchetto.

### Scostamenti dal piano, dichiarati

1. **Un test esistente cambiato** (`test-messages.R`, "print notes cover the small designs"): il dataset di 10 gruppi con seme 4 cade nell'insieme che D-90 cambia (fit senza pendenza non convergente: ora GEE, quindi senza la nota sul modello marginalizzato con 10 gruppi). Seme portato a 1 e aggiunta l'attesa `fit10$method == "mmm"`. Il piano diceva che gli altri test non cambiavano.
2. **Due errori degli script di controllo, corretti dopo il primo giro** (errori dello script, non del pacchetto; definizioni dei cancelli invariate): in `c2_study_compare.R` (e per coerenza in `c1_identity_compare.R`) l'esito atteso di `auto` richiedeva sempre una GEE stimata, mentre il piano chiede "identica a `method = "gee"`", che dove la GEE fallisce è l'errore "Neither model"; e i messaggi cli vanno a capo, quindi il testo si confronta dopo aver normalizzato gli spazi.
3. **Aggiunta non prevista**: anche l'errore "Neither model gave a usable estimate" porta il record del controllo (`convergence`), come l'errore di `method = "mmm"`.
4. **Prima esecuzione completa ripetuta**: il primo giro aveva il test del punto 1 fallito (e quindi `R CMD check` con 1 errore); i numeri riportati sono del secondo giro, con le comparazioni rifatte dopo il punto 2.

### Cosa non è cambiato

Motore del modello marginalizzato, GEE e correzioni, `gee_vcov()`, `flip_test()`, costanti dello studio (ancora quelle delle tabelle di S9, da rigenerare quando lo studio avrà le tabelle con D-90), `check_coherence = TRUE` di default. Il ramo del manoscritto e `~/R/s7d_lib` non sono stati toccati.
