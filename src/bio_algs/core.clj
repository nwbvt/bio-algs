(ns bio-algs.core
  (:require [clojure.java.io :refer :all]))

(defn write-result
  "Writes the solution in the expected format"
  [values]
  (with-open [outfile (writer "result.txt")]
    (binding [*out* outfile]
      (apply print values))))
