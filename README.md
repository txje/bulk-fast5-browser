Bulk FAST5 Browser
==================

For mild interest and low-level debugging of ONT nanopore sequencing runs

I'll add more features as I need them

Usage
-----

Serves a simple Websocket API with access to the (probably HUGE) bulk FAST5 file

    python bulk_fast5_server.py <bulk_file.fast5> [--port 8080]

For me, this works very smoothly with files multiple hundreds of GB (a full 48-hr run)

Then just open bulk_fast5_client.html in any modern (ES6-compatible) browser

The file doesn't have to be served under any sort of web server since it uses only static dependencies and the Websocket
