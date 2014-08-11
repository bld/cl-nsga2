(in-package :cl-nsga2)

(defun constraint-gte (value gte)
  (- value gte))

(defun constraint-lte (value lte)
  (- lte value))
