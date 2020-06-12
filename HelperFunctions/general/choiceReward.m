function s = choiceReward(rprob, a)

 % CHOICEREWARD
 % Based on the reward probability associated with choices, the function
 % computes if the choice will yield in a reward or not.
 %
 % INPUT: 
 %       rprob : a (1Xnumber of choices) double vector containing the choice
 %                 associated reward probabilities.
 %       a     : a choice variable. The choice coding should index the corresponding
 %               reward probability.
 %
 % OUTPUT:
 %       s     : a variable indexing the selected type of reward. 
 % =========================================================================

 s = binornd(1, rprob(a))+1;

end