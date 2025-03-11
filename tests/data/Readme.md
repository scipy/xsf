# Reference Data

## exp1

```wolfram
randomNumbers1 = RandomReal[{0.0, 1.0}, 1000]; 
randomNumbers2 = RandomReal[{0, 720}, 1000]; 
randomNumbers = Join[randomNumbers1, randomNumbers2];
results = Table[{x, -ExpIntegralEi[-x]}, {x, randomNumbers}];
formattedResults = 
  Map[StringJoin[
     StringRiffle[
      ToString[NumberForm[#, 18, ScientificNotationThreshold -> {0, 1},
          NumberMultiplier -> "*",
          NumberFormat -> (Row[{#1, "e", #3}] &)
          ]] & /@ #, ","]] &, results];
Export["exp1.csv", formattedResults, "Text"]
```
