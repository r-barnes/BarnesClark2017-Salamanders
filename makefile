salamander:
	$(MAKE) -C src/
	mv src/salamander.exe ./

clean:
	rm -f src/obj/*o
	rm -f salamander.exe src/salamander.exe src/test.exe