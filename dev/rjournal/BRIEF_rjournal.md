# Preparazione dell'articolo per The R Journal su fewclust

Aggiornato il 2026-11-03. Il manoscritto dello studio è chiuso (step S11 approvato, ramo `morris-v4`, commit `b03e36d`) e cita fewclust 0.4.2 (commit `f05c914`, ramo `package`). Il primo rilascio su CRAN è preparato come 0.4.3: fit e numeri uguali alla 0.4.2.

## Cosa serve prima di inviare l'articolo

1. **fewclust deve essere su CRAN.** La rivista considera solo pacchetti già disponibili su CRAN o Bioconductor, e rifiuta senza revisione quelli senza prove di maturità: tempo su CRAN, storia dei rilasci, sviluppo attivo. Conviene inviare la 0.4.3 a CRAN subito, scrivere l'articolo in parallelo e inviarlo dopo almeno un aggiornamento su CRAN.
2. **Repository pubblico del pacchetto.** `DESCRIPTION` punta a `github.com/lucavd/fewclust`, il repository pubblico che il PI sta creando dal ramo `package`; i due vanno tenuti allineati a ogni push. Il controllo di CRAN dà una NOTE sull'URL finché il repository non è raggiungibile.
3. **Citazione dello studio.** L'articolo deve citare il metodo, descritto altrove: serve un preprint o il riferimento "submitted".

## Requisiti della rivista (verificati il 2026-11-02)

- **Strumenti e lunghezza.** Template R Markdown di rjtools (`rjtools::create_article()`); il template Quarto è annunciato ma non ancora disponibile. Fino a 20 pagine; titolo e abstract in testo semplice, abstract di non più di 250 parole; output prodotti dal codice.
- **Struttura per un pacchetto:**
  - introduzione con il problema, i pacchetti concorrenti e il vuoto colmato;
  - metodo solo quanto serve;
  - scelte di progetto;
  - flusso tipico;
  - un caso d'uso;
  - dettagli tecnici;
  - discussione.
- **Revisione.** Servono vantaggi dimostrati rispetto ai concorrenti, dipendenze contenute e il controllo di versione collegato nel `DESCRIPTION`.
- **Controlli prima dell'invio:** `rjtools::initial_check_article()` (pacchetto su CRAN, pacchetti citati disponibili, ortografia, DOI della bibliografia).

## Concorrenti da confrontare (stato su CRAN verificato il 2026-11-02)

- **Su CRAN:**
  - geepack (motore GEE usato da fewclust);
  - geesmv 1.3;
  - geessbin 1.0.2;
  - glmtoolbox 0.1.12;
  - geeCRT;
  - clubSandwich 0.7.0;
  - GLMMadaptive 0.9-7 (coefficienti marginali dai modelli misti);
  - flipscores 1.3.2 (test sign-flip degli autori del metodo).
- **Non su CRAN:** MMLB (solo GitHub).

Confronti numerici da fare nell'articolo:
- `gee_vcov()` contro geesmv e clubSandwich. Il supplemento del paper documenta che geesmv applica a Kauermann-Carroll una radice quadrata simmetrica a una matrice non simmetrica, e che con correlazione exchangeable geesmv e fewclust coincidono per Mancl-DeRouen e Fay-Graubard.
- Il modello marginalizzato contro MMLB: riferimenti già nei test.
- `flip_test()` contro flipscores.
- I tempi.

## Struttura proposta (15-18 pagine)

1. Introduction.
2. Models and inference in brief:
   - modello marginalizzato con controllo di convergenza (regola del fit di riserva solo se converge);
   - t con N - p gradi di libertà;
   - correzioni della GEE;
   - test sign-flip;
   - criterio di stima usabile.
3. Design of the package:
   - una funzione con default dallo studio;
   - condizioni con classe;
   - interfaccia tidyverse;
   - perché i motori sono nel pacchetto.
4. A typical workflow.
5. Use case: dati pubblici (per esempio il trial respiratorio di geepack, già usato nel paper) o dati simulati.
6. Comparison with existing software.
7. What the simulation study found: numeri letti dalle costanti del pacchetto.
8. Implementation details.
9. Discussion and limitations.

## Materiale già pronto

- Output dimostrativi della 0.4.0 (stampa, `summary()`, `tidy()`, `glance()`, GEE, `flip_test()`, tempi) negli artifact della sessione del 2026-11-02. Vanno rigenerati con la versione su CRAN.
- Collezione Zotero "MMM" (84 titoli distinti). Mancano Heagerty 1999 (`10.1111/j.0006-341X.1999.00688.x`), Hemerik, Goeman e Finos 2020 (`10.1111/rssb.12369`) e i riferimenti dei pacchetti concorrenti.

## Regole per la sessione

- Articolo in una cartella `rjournal/` del nuovo repository, o alla radice del ramo `package` (fuori dal pacchetto costruito).
- `~/MMM` resta in sola lettura e non si installa nulla in `~/R/s7d_lib`.
- Sul server mancano rjtools, distill, bookdown, kableExtra, flipscores, geessbin, glmtoolbox e GLMMadaptive: vanno installati nella libreria degli strumenti `~/MMM-package/.Rlib`.
- Nei job `pdflatex` va aggiunto al `PATH` (`$HOME/bin`).
