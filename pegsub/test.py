# def application(env, start_response):
#     start_response('200 OK', [('Content-Type','text/html')])
#     return ["Hello World"] # python3
def application(env, start_response):
    start_response('200 OK', [('Content-Type','text/html')])
    return [b"Hello World"] # python3
