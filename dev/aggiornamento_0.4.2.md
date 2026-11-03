# fewclust 0.4.2: costanti dello studio dalle tabelle definitive di S11

Preparato il 2026-11-03 su richiesta dell'agente principale: le tabelle dello step S11 sono definitive (pipeline di chiusura `final` e passaggio del testo `final2` sulla copia di lavoro di `~/MMM`, ramo `morris-v4`, non ancora committate; registro dei numeri `tables/V5_S7f_numbers.csv`, 3.376 numeri, 0 discordanze). Le costanti del pacchetto che citano risultati della simulazione (`R/constants.R`, lette dalla vignetta "What the simulation study found") vengono ancora dalle tabelle di S9 e si rigenerano da quelle di S11. Si parte dalla 0.4.1 (commit `45c7946`). `~/MMM` si legge soltanto.

**Stato di questo file.** Le sezioni fino a "Controlli" sono scritte e committate prima di qualsiasi modifica. Gli esiti si aggiungono alla fine.

## Cosa cambia

1. `dev/make_constants.R` legge gli stessi quattro file di `~/MMM/tables` (`V5_S7f_table4.csv`, `V5_S7f_rule_dp.csv`, `V5_S7f_cp_S7d_dp_rates.csv`, `V5_S7f_h2h_boot.csv`, ora nella versione di S11, con md5 uguale a quello del registro dove il registro li cita) e il registro `V5_S7f_numbers.csv`. Ogni valore di `R/constants.R` porta accanto, in commento, il file sorgente con la riga (strategia, regione o punto del disegno) e l'id del registro che cita quella cella (stesso file, stesso filtro, stesso campo), oppure "no registry id" se il registro non la cita.
2. Nessun numero scritto a mano: anche i valori finora digitati nel generatore (10, 15 e 145 gruppi, 10 occasioni, 4 meccanismi, la soglia di 3 occasioni e 25 gruppi della regione "piccola") si leggono da `V5_S7f_rule_dp.csv` (colonne `N`, `ntime`, `dgp`, `small`).
3. Vignetta: le frasi che oggi affermano un risultato a parole (il modello marginalizzato "mai" sopra 7,5%; le altre correzioni sopra 7,5% "più spesso"; il test sign-flip "ha mantenuto il livello"; il margine di 1 punto nella didascalia) si calcolano dalle costanti, o si tolgono se non sono più vere.
4. Versione: 0.4.2 se cambia qualcosa che l'utente vede (un numero della vignetta), con una voce in NEWS; altrimenti resta 0.4.1 con un commit di documentazione. Atteso 0.4.2: con D-90 la quota di risultati usabili del modello marginalizzato scende (da 90,9% a 84,0% secondo D-90).
5. Nient'altro cambia: codice in `R/` salvo `R/constants.R`, test, `fit_fewclust()`, messaggi.

## Controlli (cancelli)

1. **Costanti = tabelle** (`dev/checks_0.4.2/c1_constants.R`): ogni valore di `.fc_study` della libreria installata è uguale alla cella del file sorgente (interi e testi `identical()`, reali `all.equal(tolerance = 1e-12)`) e, dove c'è un id, al valore del registro; md5 dei file sorgente uguale a quello scritto nell'intestazione e a quello del registro; struttura (nomi, colonne, righe) uguale alla 0.4.1, così la vignetta si costruisce.
2. **Solo costanti e documentazione**: `git diff` rispetto alla 0.4.1 tocca in `R/` solo `R/constants.R`.
3. **Test, lint, `R CMD check --as-cran`** come per la 0.4.1: 0 test falliti, 0 lint in `R/`, 0 errori e 0 warning.
4. **Elenco delle costanti cambiate** (vecchio -> nuovo, id del registro) in `dev/checks_0.4.2/out/constants_diff.csv`, letto dalle due versioni di `R/constants.R`.

## Esiti

Eseguiti il 2026-11-03 con `bash dev/checks_0.4.2/run_checks.sh <R/constants.R della 0.4.1>` sul server (R 4.6.1, un processo R); uscite in `dev/checks_0.4.2/out/`. Sorgenti lette dalla copia di lavoro di `~/MMM/tables` (non committata): md5 scritti nell'intestazione di `R/constants.R`, uguali a quelli del registro `V5_S7f_numbers.csv` per `V5_S7f_table4.csv` e `V5_S7f_h2h_boot.csv` (il registro non cita `V5_S7f_rule_dp.csv` e `V5_S7f_cp_S7d_dp_rates.csv`).

| Cancello | Esito | Dettaglio |
|---|---|---|
| 1. Costanti = tabelle | PASS | 375 valori (table4 285, ar1_fg_hi 18, h2h 63, scalari 9) uguali alle celle delle tabelle; 159 valori portano almeno un id del registro (168 id), tutti uguali al valore del registro; gli altri 216 portano "no registry id" (colonne che il registro non cita: `n_dp`, `n_main`, ..., `margin`, i punti AR(1) con Fay-Graubard, gli scalari); md5 sorgenti = intestazione = registro; struttura uguale alla 0.4.1 |
| 2. Solo costanti e documentazione | PASS | in `R/` cambia solo `R/constants.R` |
| 3. Test, lint, `R CMD check --as-cran` | PASS | 13 file, 569 attese, 0 fallite, 0 errori, 0 avvisi; 0 lint in `R/`; 0 errori, 0 warning, 1 NOTE (HTML Tidy non installato) |
| 4. Elenco delle costanti cambiate | fatto | `out/constants_diff.csv`: 25 valori, tutti del modello marginalizzato (table4, riga `mmm_cc`: `lo_all`, `lo_N10`, `pu`, `upow`, `pu_small`; h2h: `mad_mmm` e `diff` nelle 9 celle, `class` in 2) |

Versione 0.4.2: i valori cambiati compaiono nella vignetta. Le frasi della vignetta che dipendono dai numeri (modello marginalizzato mai sopra 7,5%; le altre GEE sopra 7,5% più spesso di AR(1) con Fay-Graubard; il test sign-flip mantiene il livello; il margine della didascalia) ora si calcolano dalle costanti; con le tabelle di S11 restano tutte vere.

### Scostamenti dal piano, dichiarati

1. Nella prima prova di `c1_constants.R` (prima del giro registrato) la lettura dei valori testuali con spazi ("about the margin") era sbagliata nello script di controllo; corretta, il generatore verifica comunque ogni id contro il registro.
2. La voce NEWS descrive i cambiamenti senza numeri (rimanda alla vignetta e a `constants_diff.csv`), per non scrivere numeri a mano.
3. La frase della vignetta "no method with an interval was both available and calibrated" (2-3 occasioni, al massimo 25 gruppi) resta testo: non si calcola dalle costanti salvate. Con D-90 il modello marginalizzato è disponibile meno spesso, quindi la frase non diventa falsa.
4. `R CMD check` usa lo script della 0.4.1 (`dev/checks_0.4.1/c7_rcmdcheck.sh`), che scrive "GATE D" nel riassunto.
