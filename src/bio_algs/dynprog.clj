(ns bio-algs.dynprog)

(defn find-change
  "given a money value and a list of coins, return the minimum number of coins needed to make that amount"
  [value coins]
  (if (zero? value) 0
    (inc (apply min (for [coin coins :when (<= coin value)]
                      (find-change (- value coin) coins))))))
