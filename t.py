from ipaddress import ip_network

net = ip_network("192.192.192.192/255.255.254.0",0)
print(net)
print(net.num_addresses)
print(net.network_address)
enumerate(net,)
con

