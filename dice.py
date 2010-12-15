#!/usr/bin/env python

'''
Sundry dice simulations, with brute force and Monte Carlo analysis.

URL:            <http://www.sailmaker.co.uk/newfiles/dice.py>

Maintainer:     Richard Buckle <mailto:richardb@sailmaker.co.uk>

Licence:        Public domain

Usage:          See main()

Compatibility:  Python 2.4 or later
'''


###### imports ######
import random
import math


###### constants ######
__author__      = "Richard Buckle <mailto:richardb@sailmaker.co.uk>"
__status__      = "stable"
__version__     = "1.0.3"
__revision__    = "1.0.3"
__date__        = "2009.01.16"
__copyright__   = "Public domain"


##### class Histogram #####
class Histogram(dict):
    '''
    Subclass of dict to represent a histogram.
    
    The dict.key is the bucket, e.g. the value obtained by rolling some dice
    under whatever set of rules.
    The corresponding dict.value is the observed frequency of that bucket.
    
    To add a sample, use something like: 
        histogram[bucket] = histogram[bucket] + 1
    '''
    
    def __init__(self, minval, maxval):
        'Create all buckets and initialize them to 0.'
        for i in xrange(minval, maxval + 1):
            self[i] = 0
            
    def numSamples(self):
        'Return the number of samples.'
        return sum(self.itervalues())
        
    def sigma(self):
        'Return the weighted sum of values.'
        return sum([val * freq for val, freq in self.iteritems()])
        
    def sigmaSquared(self):
        'Return the weighted sum of squared values.'
        return sum([val * val * freq for val, freq in self.iteritems()])
        
    def frequency(self, k):
        'Return the frequency of key k as a percentage.'
        return 100.0 * self.get(k, 0) / self.numSamples()
        
    def mean(self):
        'Return the arithmetic mean.'
        return 1.0 * self.sigma() / self.numSamples()
        
    def popStdDev(self):
        'Return the population standard deviation.'
        n = self.numSamples()
        sigma = self.sigma()
        sigmaSq = self.sigmaSquared()
        variance = 1.0 * (n * sigmaSq - sigma * sigma) / (n * n)
        return math.sqrt(variance)
        
    def sampleStdDev(self):
        'Return the sample standard deviation.'
        n = self.numSamples()
        sigma = self.sigma()
        sigmaSq = self.sigmaSquared()
        variance = 1.0 * (n * sigmaSq - sigma * sigma) / (n * (n-1))
        return math.sqrt(variance)
        
    def printTabbed(self):
        'Print tab-delimited, for import by apps such as Excel.'
        for bucket in sorted(self.keys()):
            print "%s\t%s" % (bucket, self[bucket])
        
    def dump(self):
        'Print results and summary statistics.'
        self.printTabbed()
        print 'mean = %.2f' % self.mean()
        print 'sampleStdDev = \t%.2f' % self.sampleStdDev()
        print 'popStdDev = \t%.2f' % self.popStdDev()
        print


##### class MultiHistogramTabulator #####
class MultiHistogramTabulator(dict):
    '''
    A class to tabulate multiple histograms together for output,
    with headers to distinguish each one.
    
    Tabulates the union of the buckets of each component histogram,
    together with the summary statistics of each.
    '''
    
    def __init__(self, histograms, headers):
        'Construct from a list of histograms and a list of headers.'
        self.headers = headers[:]
        self.means = []
        self.sampleStdDevs = []
        self.popStdDevs = []
        
        # accumulate set of buckets and lists of summary statistics
        allBuckets = set()
        for h in histograms:
            # accumulate buckets
            allBuckets.update(h.iterkeys())
            # accumulate summary statistics
            self.means.append(h.mean())
            self.sampleStdDevs.append(h.sampleStdDev())
            self.popStdDevs.append(h.popStdDev())
        
        # accumulate percentage frequencies over all buckets
        for bucket in allBuckets:
            self[bucket] = [h.frequency(bucket) for h in histograms]
            
    def printTabbed(self):
        'Print tab-delimited, for import by apps such as Excel.'
        # print headers
        print '\t%s' % '\t'.join([str(x) for x in self.headers])
        # print bucket value and list of frequencies
        for val, freqs in self.items():
            print '%s\t%s' % ( val, '\t'.join(['%.1f' % x for x in freqs]) )
        # finish with a blank line
        print
        
    def dump(self):
        'Print space-delimited, for command line etc.'
        # print headers
        print '\t%s' % '\t'.join([str(x).center(5) for x in self.headers])
        # print bucket value and list of frequencies
        for val, freqs in self.items():
            print '%s\t%s' % ( val, '\t'.join(['%5d' % x for x in freqs]) )
        # print lists of summary statistics
        print 'avg\t%s' % '\t'.join(['%5.2f' % x for x in self.means])
        print 'ssd\t%s' % '\t'.join(['%5.2f' % x for x in self.sampleStdDevs])
        print 'psd\t%s' % '\t'.join(['%5.2f' % x for x in self.popStdDevs])
        # finish with a blank line
        print
        
      
##### class Dice #####
class Dice(object):
    '''
    Base class for rolling dice; sums the pips on the dice.
    
    To implement a different rule set, it is usually only necessary to override
    roll(), for Monte Carlo analysis, and fromState(), for brute force analysis.
    
    Subclasses that implement "roll again" rule sets (e.g. d20 critical hits and 
    the StoryTelling system) will usually need to allocate additional "shadow dice" 
    in self.state for the potential second roll in order to do a full brute force analysis. 
    
    Be aware that this will exponentially increase the time taken by brute force analysis.
    In such cases, a Monte Carlo approach is perhaps the better option.
    '''

    def __init__(self, num, sides, multiplier=1):
        'Initialise from number of dice, number of sides and overall multiplier.'
        if num < 1:
            raise ValueError('No dice given!')
        self.num = num
        self.sides = sides
        self.multiplier = multiplier
        self.minVal = num * multiplier 
        self.maxVal = num * multiplier * sides 
        self.state = [1] * num # list of dice values, for brute force analysis
        
    def __iter__(self):
        'Return an iterator, for use by bruteforce().'
        return self
        
    def roll(self):
        '''
        Evaluate one random roll of these dice.
        
        Subclasses will usually need to override this.
        '''
        def accumulateRoll(oldVal, dummy, sides=self.sides):
            # Roll another of these dice and add it to the result.
            return oldVal + random.randint(1, sides)
        return reduce(accumulateRoll, range(self.num+1)) * self.multiplier
        
    def _normaliseState(self):
        '''
        Used internally by next().
        
        If any die in self.state has been advanced too far, 
        turn it back to 1 and turn the next higher die, recursively.
        
        Raises StopIteration when the final state is reached, which is indicated
        by the leftmost die reaching its maximum value.
        '''
        # examine dice from right to left
        for i in xrange(len(self.state)-1, -1, -1):
            if self.state[i] > self.sides:
                # this die has been turned too far
                if i == 0:
                    # leftmost die turned too far, so we're done
                    raise StopIteration
                # turn this die back to 1
                self.state[i] = 1 
                # turn the next higher die
                self.state[i-1] += 1 
                # re-examine
                self._normaliseState()
    
    def fromState(self):
        '''
        Evaluate the current brute force state:
        the base class implementation is to sum pips on the dice.
        
        Subclasses will usually need to override this.
        '''
        return sum(self.state)

    def next(self):
        '''
        Used internally by bruteforce().
        
        Return our evaluation of the current brute force state, also advancing it.
        '''
        self._normaliseState() # can raise StopIteration
        # compute current state
        result = self.fromState()
        # turn the rightmost die, ready for the next call
        self.state[-1] += 1
        return result

    def montecarlo(self, trials):
        'Perform a monte carlo simulation and return the histogram.'
        histogram = Histogram(self.minVal, self.maxVal)
        for dummy in xrange(trials):
            # random roll
            roll = self.roll()
            histogram[roll] = histogram[roll] + 1
        return histogram
        
    def bruteforce(self):
        'Perform a brute force enumeration and return the histogram.'
        histogram = Histogram(self.minVal, self.maxVal)
        for roll in self:
            histogram[roll] = histogram[roll] + 1
        return histogram
    
    
##### class DiceDiscard #####
class DiceDiscard(Dice):
    '''
    Like Dice, but a given number of lowest values are discarded.
    
    Sometimes used for character generation in d20 systems.
    '''
    
    def __init__(self, num, sides, multiplier=1, discards=0):
        '''
        Initialise from number of dice, number of sides, overall multiplier 
        and number of discardable dice.
        
        Requires discards <= num.
        '''
        Dice.__init__(self, num, sides, multiplier)
        self.discards = discards
        self.minVal = (num - discards) * multiplier 
        self.maxVal = (num - discards) * multiplier * sides 
        if self.discards > num:
            raise ValueError('More discards than dice.')
        
    def _evaluateWithDiscards(self, dice):
        'Evaulate a list of dice by discarding the lowest values.'
        sortedDice = dice[:]
        sortedDice.sort()
        dicePicked = sortedDice[self.discards:]
        return sum(dicePicked) * self.multiplier
        
    def roll(self):
        'Evaluate one random roll of these dice.'
        dice = []
        for dummy in xrange(self.num):
            dice.append(random.randint(1, self.sides))
        return self._evaluateWithDiscards(dice)
        
    def fromState(self):
        'Evaluate the current brute force state.'
        return self._evaluateWithDiscards(self.state)
    
    
##### class StoryTellerDice #####
class StoryTellerDice(Dice):
    '''
    Dice rolls for systems such as White Wolf's "StoryTelling" system.
    <http://en.wikipedia.org/wiki/Storyteller_System>
    
    Like Dice, but result is number of dice with pips >= threshold.
    Dice with pips >= rollAgain can be rolled again,
    for +1 to the result if the pips on the reroll >= threshold.
    
    TODO: this implementation only allows roll-again only once, whereas
    White Wolf's "StoryTelling" system allows unlimited roll-again.
    '''
    
    def __init__(self, num, sides=10, threshold=8, rollAgain=10):
        '''
        Initialise from number of dice , number of sides, 
        threshold for success and threshold to roll again.
        
        TODO: Roll again is once-only, contrary to White Wolf's system.
        
        Requires rollAgain >= threshold.
        '''
        #sanity check
        if rollAgain < threshold:
            raise ValueError('rollAgain is less than threshold.')
        
        Dice.__init__(self, num, sides, 1)
        self.threshold = threshold
        self.rollAgain = rollAgain
        self.minVal = 0 
        self.maxVal = num * 2 # all roll agains succeed
        
        # Allocate "shadow dice" for "roll agains" in brute force analysis:
        # the initial rolls are self.state[0:self.num]
        # the "roll again" roll for any self.state[i] is self.state[i + self.num].
        #
        # This is obviously a very dumb way to do brute force analysis
        # given that a roll has only three outcomes, rather than self.sides outcomes, 
        # but for this cut I'll focus on correctness and maybe develop a 
        # smarter algorithm later.
        # 
        # Preliminary analysis is that a Monte Carlo of 50000 runs didn't 
        # significantly differ from brute force for a dice pool of 10 dice.
        self.state = [1] * (num * 2)
        
    def roll(self):
        'Evaluate one random roll of these dice.'
        successes = 0
        for dummy in xrange(self.num):
            diceRoll = random.randint(1, self.sides)
            if diceRoll >= self.threshold:
                # success
                successes += 1
                if diceRoll >= self.rollAgain:
                    # roll again
                    if random.randint(1, self.sides) >= self.threshold:
                        successes += 1
        return successes
        
    def fromState(self):
        '''
        Evaluate the current brute force state.
        rollAgain giving another >= self.rollAgain is ignored.
        '''
        successes = 0
        for i in xrange(0, self.num):
            diceRoll = self.state[i]
            if diceRoll >= self.threshold:
                # success
                successes += 1
                if diceRoll >= self.rollAgain:
                    # roll again, shadow die is self.state[i + self.num]
                    if self.state[i + self.num] >= self.threshold:
                        successes += 1
        return successes
    
    
##### class RiskDice #####
class RiskDice(Dice):
    '''
    Dice rolls for the game of Risk.
    
    Attacker rolls attackDice d6, defender rolls defendDice d6,
    where defendDice <= attackDice.
    Each player's dice are sorted and compared.
    If the attacker's die beats the defender's die, the defender loses,
    otherwise the attacker loses.
    '''
    
    def __init__(self, attackDice, defendDice):
        '''
        Initialise from number of attacking and defending dice.
        
        Requires defendDice <= attackDice.
        '''
        Dice.__init__(self, 1, 6)
        self.attackDice = attackDice
        self.defendDice = defendDice
        if defendDice > attackDice:
            raise ValueError('More defend dice than attack dice')
        
    def sortedRoll(self, num):
        'Return one random roll of num dice, sorted descending.'
        sortedDice = []
        for dummy in xrange(num):
            sortedDice.append(random.randint(1, self.sides))
        sortedDice.sort()
        sortedDice.reverse()
        return sortedDice

    def roll(self):
        '''
        Return the outcome of one random battle with these dice
        as a tuple (numDefenderWins, numAttackerWins)
        '''
        attackRoll = self.sortedRoll(self.attackDice)
        defendRoll = self.sortedRoll(self.defendDice)
        
        def fightFn(attack, defend):
            '''
            Return the outcome of one random comparison of these dice
            as a tuple (defenderWins, attackerWins)
            '''
            if attack is None or defend is None:
                return None
            elif attack > defend:
                # attacker wins
                return (0, 1)
            else:
                # defender wins
                return (1, 0)
        
        results = map(fightFn, attackRoll, defendRoll)
        
        def summarize(pair1, pair2):
            'add two pairs together, handling None'
            if pair1 is None: 
                return pair2
            if pair2 is None: 
                return pair1
            return (pair1[0] + pair2[0], pair1[1] + pair2[1])
        
        results = reduce(summarize, results)
        return results
        
    class Results(object):
        '''
        Square zero-based Matrix of [defender wins][attacker wins] 
        '''
        def __init__(self, defendDice):
            'Create matrix and initialize to all zeros'
            self.matrix = []
            self.numResults = 0
            row = [0] * (defendDice + 1)
            for i in xrange(defendDice + 1):
                self.matrix.append(row[:])
                
        def addResult(self, defenderWins, attackerWins):
            'Add the result to the matrix'
            self.matrix[defenderWins][attackerWins] += 1
            self.numResults += 1
            
        def _dims(self):
            'xrange for our (square) dimension'
            return xrange(0, len(self.matrix))

        def dump(self, percentages=True):
            '''
            Print results in a matrix of the form
                              Attacker wins
            Defender wins   +-----0-----1-----2
                          0 |   0.0   0.0  37.6
                          1 |   0.0  33.1   0.0
                          2 |  29.3   0.0   0.0

            
            If percentages=True, display percentages,
            otherwise display the number of rolls.
            '''
            # print header
            print '                  Attacker wins'
            print 'Defender wins   +-----' + ('-----'.join(['%u' % i for i in self._dims()]))
            # print table rows
            for defenderWins in self._dims():
                row = self.matrix[defenderWins]
                if percentages:
                    rowText = ''.join([('%6.1f' % (100.0 * row[i] / self.numResults)) for i in self._dims()])
                else:
                    rowText = ''.join(['%6d' % row[i] for i in self._dims()])
                print '              %u |%s' % (defenderWins, rowText)
            print

    def montecarlo(self, trials):
        '''
        Perform a monte carlo simulation and return the histogram
        '''
        results = self.Results(self.defendDice)
        for i in xrange(trials):
            roll = self.roll()
            results.addResult(roll[0], roll[1])
        return results
    
    
##### class ArcanumDice #####
class ArcanumDice(Dice):
    '''
    Dice rolls for the game found in Arcanum.
    
    The player rolls 2d6.
    On the first roll, 7 is a win, 2 is a loss and any other score is remembered as the 'mark',
    in which case the player rolls again.
    On subsequent rolls, the mark is a win, 2 and 7 are both losses,
    and anything else means roll again.
    '''
    
    def __init__(self):
        '''
        Initialise from the fixed number of dice specified by the rules.
        
        The histogram is +1 for win, -1 for loss and (internal only) 0 for roll again.
        '''
        Dice.__init__(self, 2, 6)
        self.mark = 0
        self.minVal = -1
        self.maxVal = 1
        
    def _evaluate(self):
        'Evaluate one roll in the match.'
        pips = Dice.roll(self)
        if self.mark == 0: # first roll 
            if pips == 7:
                return 1 # win
            if pips == 2:
                return -1 # loss
            self.mark = pips
            return 0 # roll again
        else:
            if pips == 7 or pips == 2:
                return -1 # loss
            if pips == self.mark:
                return 1 # win
            return 0 # roll again
    
    def roll(self):
        'Evaluate one match.'
        self.mark = 0
        result = 0
        while result == 0:
            result = self._evaluate()
        return result
        
    def bruteforce(self):
        print 'Brute force is not supported for the Arcanum dice rules.'

    
##### main #####
def main():
    '''Examples of usage'''
    
    if False:
        print '4d6 montecarlo'
        dice = Dice(4, 6)
        histogram = dice.montecarlo(10000)
        histogram.dump()
        
        print '2d6 * 2 montecarlo'
        dice = Dice(2, 6, 2)
        histogram = dice.montecarlo(10000)
        histogram.dump()
    
    if False:
        print '3d6 bruteforce'
        dice = Dice(3, 6)
        histogram = dice.bruteforce()
        histogram.dump()

    if False:
        print '3d6 montecarlo'
        dice = Dice(3, 6)
        histogram = dice.montecarlo(10000)
        histogram.dump()

    if False:
        print '4d6 discard 1 montecarlo'
        dice = DiceDiscard(4, 6, 1, 1)
        histogram = dice.montecarlo(10000)
        histogram.dump()

    if False:
        print '4d6 discard 1 bruteforce'
        dice = DiceDiscard(4, 6, 1, 1)
        histogram = dice.bruteforce()
        histogram.dump()

    if False:
        print 'Risk montecarlo'
        for armies in ( (3, 2), (3, 1), (2, 2), (2, 1), (1, 1) ):
            print '%d on %d' % armies
            dice = RiskDice(*armies)
            results = dice.montecarlo(10000)
            results.dump()
    
    if False:
        print 'Storyteller montecarlo'
        histograms = []
        maxDice = 5
        for i in xrange(1, maxDice + 1):
            dice = StoryTellerDice(i)
            histogram = dice.montecarlo(50000)
            histograms.append(histogram)
        multihisto = MultiHistogramTabulator(histograms, range(1, maxDice + 1))
        multihisto.dump()
    
    if False:
        print 'Storyteller bruteforce'
        histograms = []
        maxDice = 2
        for i in xrange(1, maxDice + 1):
            dice = StoryTellerDice(i)
            histogram = dice.bruteforce()
            histograms.append(histogram)
        multihisto = MultiHistogramTabulator(histograms, range(1, maxDice + 1))
        multihisto.dump()

    if True:
        print 'Arcanum montecarlo'
        dice = ArcanumDice()
        histogram = dice.montecarlo(5)
        histogram.dump()


if __name__ == "__main__":
    main()
