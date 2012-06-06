;;; sand-pile.lisp --- 2D CA simulating a sand pile

;; Copyright 2012 Eric Schulte

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; This program implements the 2D CA used to investigate
;; self-organized criticality through the simulation of a simplified
;; sand-pile model in the following publication.

;; @article{bak1988self,
;;   author={Bak, P. and Tang, C. and Wiesenfeld, K. and others},
;;   title={Self-organized criticality},
;;   journal={Physical review A},
;;   year={1988},
;;   volume={38},
;;   number={1},
;;   pages={364--374}
;; }

;;; Code:
(require :alexandria) (require :metabang-bind)
(defpackage #:sand-pile (:use :common-lisp :alexandria :metabang-bind))
(in-package :sand-pile)

(defvar *size* 20)
(defvar *dim* 2)
(defvar *grid* (make-array (make-list *dim* :initial-element *size*)
                           :element-type :number))

(defun settle (rule)
  "Apply RULE to all pairs of neighbors in the `*grid*'.
RULE should take two neighbor values at time step n as input and
return their two values at time step n+1.")

(defun init (rule)
  "Use RULE to initialize the `*grid*'.
RULE should take grid coordinates and return the starting value."
  )

(defun run (times rule &aux (history (list *grid*)))
  "Run the system for TIMES time steps using update RULE and `settle'.
Return a list of the states of the system at each time step."
  (dotimes (n times history) (settle rule) (push *grid* history)))


;;; Utility
(defun range (n)
  (loop for i upto (1- n) collect i))

(defun cross (dimensions)
  (if (cdr dimensions)
      (let ((this (range (car dimensions))))
        (mapcan (lambda (rst) (mapcar (lambda (el) (cons el rst)) this))
                (cross (cdr dimensions))))
      (mapcar #'list (range (car dimensions)))))

(defun neighbors (dim)
  (remove-if (lambda (dim) (some (lambda (n) (or (< n 0) (< *size* n))) dim))
             (mapcar (lambda (off) (map 'list #'+ off dim))
                     (mapcar (lambda (off) (mapcar #'1- off))
                             (remove (make-list *dim* :initial-element 1)
                               (cross (make-list *dim* :initial-element 3))
                               :test #'tree-equal)))))

;;; sand-pile.lisp ends here
