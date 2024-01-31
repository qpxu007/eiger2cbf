import fabio

files = ['insu6_1_000080.cbf', 'ref/insu6_1_000080.cbf']
f1 = fabio.open(files[0], "r")
f2 = fabio.open(files[1], "r")

d= f1.data-f2.data
print(d.max(), d.min())

cbf_i = fabio.cbfimage.CbfImage(f1.data-f2.data)
cbf_i.save("x_001.cbf")
