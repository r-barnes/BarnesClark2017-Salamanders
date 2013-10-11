salamander:
	$(MAKE) -C src/
	mv src/salamander.exe ./

clean:
	rm -f src/obj/*o
