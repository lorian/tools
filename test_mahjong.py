# Unit tests for mahjong.py
import mahjong
from mahjong import tile
from mahjong import hand

def test_tile():
	test_tile = tile('suit','number')
	# __init__()
	assert test_tile.suit == 'suit'
	assert test_tile.number == 'number'
	# __repr__()
	assert eval(repr(test_tile)) == test_tile
	# __eq__()
	assert test_tile == tile('suit','number')
	assert test_tile != tile('suit','number2')
	assert test_tile != tile('suit2','number')

def test_hand():
	# __init__(), draw()
	faces = {'wind':['N','S','E','W'],'dragon':['R','G','W']}
	deck = [tile(S,N) for S in faces for N in faces[S]]
	deck.append(tile('flower',0)) # [tile('wind','N'), tile('wind','S'), tile('wind','E'), tile('wind','W'), tile('dragon','R'), tile('dragon','G'), tile('dragon','W'), tile('flower','0')]
	test_hand = hand(5,deck)
	assert test_hand.public == [tile('flower',0)]
	assert test_hand.private == [tile('dragon','W'), tile('dragon','G'), tile('dragon','R'), tile('wind','W'), tile('wind','E')]
	# discard()
	assert tile('dragon','G') in test_hand.private
	test_hand.discard(tile('dragon','G'))
	assert tile('dragon','G') not in test_hand.private
	# place()
	test_hand.place(tile('dragon','W'), tile('wind','E'), tile('dragon','R'), tile('wind','W'))
	assert test_hand.public == [tile('flower',0),tile('dragon','W'), tile('wind','E'), tile('dragon','R'), tile('wind','W')]
	assert test_hand.private == [tile('wind','S')]


# convert mahjong into __main__
