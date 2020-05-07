for ((i = 16; i <= 4096 ; i*=2))
do
	ROWA=$((i))
	COLA=$((i))
	COLB=$((i))
	./OutputProg $ROWA $COLA $COLB
done