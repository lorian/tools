# Plays mahjong using objects

import itertools
import random
import collections

class tile:
	'''
	Tile object contains suit and number of tile.
	Suit can be: bamboo, character, dot, wind, dragon, or flower.
	Wind's number is N,S,E or W. Dragon's number is R,G, or W. Flower always has a number of 0.
	'''
	def __init__(self, suit, number):
		self.suit = suit
		self.number = number

	def __eq__(self,other):
		return self.suit == other.suit and self.number == other.number

	def __hash__(self):
		return hash((self.suit,self.number))

	def __repr__(self):
		# assert eval(repr(testtile)) == testtile
		return "tile('{0}','{1}')".format(self.suit,self.number)

	def __str__(self):
		 return '{0} {1}'.format(self.number,self.suit)

#hand object, can report public tiles
class hand:
	'''
	Hand object contains reference to deck, and however many tiles are declared at initiation.
	Tiles are divided into public and private.
	'''
	def __init__(self,size,deck):
		self.deck = deck
		self.private = []
		self.public = [] #list of tuples or single flowers
		# Draw starting hand
		for n in range(size):
			self.draw()
		print self

	def __str__(self):
		return "Hand: {0}\nPlayed: {1}".format(self.private,self.public)

	def draw(self):
		'''Pops new tile from deck. If flower, puts tile in public and replaces. Returns new tile after adding to hand.'''
		new_tile = self.deck.pop()
		if new_tile.suit == 'flower':
			print "Playing flower."
			self.public.append(new_tile)
			self.draw()
		else:
			self.private.append(new_tile)
		return new_tile

	def discard(self,tile):
		'''Discards tile from hand, returns discarded tile.'''
		self.private.remove(tile)
		return tile

	def place(self,t1,t2,t3,t4=None):
		'''Play 3 or 4 tile set in public; replace 4th tile.'''
		self.public.extend([t1,t2,t3])
		self.private.remove(t1)
		self.private.remove(t2)
		self.private.remove(t3)
		if t4:
			self.public.append(t4)
			self.private.remove(t4)
			return self.draw() # replace 4th tile
		return

class player:
	'''
	Has single hand object.
	'''
	def __init__(self,fresh_hand):
		self.myhand = fresh_hand

	def turn(self):
		'''Player takes standard turn'''
		new_tile = self.myhand.draw()
		print "Drew: {0}".format(new_tile)
		self.myhand.discard(new_tile)
		print "Discarded: {0}".format(new_tile)

	def find_sets(self):
		'''Find pairs, trios, and kongs'''
		print self.myhand.private
		match = collections.Counter(self.myhand.private)
		print match

	def find_runs(self):
		'''Find sequential runs of two or three tiles'''

#player object, contains playing AI
#board/game object

# Create shuffled deck with 4 copies of each tile, so each draw can be simulated with deck.pop()
suits = ['bamboo','character','dot']
faces = {'wind':['N','S','E','W'],'dragon':['R','G','W'],'flower':[0,0]}
deck = [tile(S,N) for N in range(1,10) for S in suits]*4 + [tile(S,N) for S in faces for N in faces[S]]*4
random.shuffle(deck)
P1 = player(hand(16,deck))
P2 = player(hand(16,deck))
P1.find_sets()

