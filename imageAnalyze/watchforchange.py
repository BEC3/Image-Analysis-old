import time
import sys, os
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
# import AnalyzeAIA


class MyHandler(PatternMatchingEventHandler):
	def __init__(self):
		super( MyHandler, self ).__init__()
		# self.f = functionToCall
		# self.parent = parent

	patterns = ["*.aia"]

	def process(self, event):
		"""
		event.event_type
			'modified' | 'created' | 'moved' | 'deleted' 
		event.is_directory
			True | False
		event.src_path
			path/to/observed/file

		"""
		# the file will be processed there
		print event.src_path, event.event_type  # print now only for debug
		# self.f(self.parent)
		# self.showImg()
		# os.system('python AnalyzeAIA.py')
	def on_modified(self, event):
		self.process(event)
	def on_created(self, event):
		self.process(event)

if __name__ == '__main__':
    args = sys.argv[1:]
    observer = Observer()
    observer.schedule(MyHandler(), path = args[0] if args else '.')
    observer.start()

    try:
	while True:
	   time.sleep(1)
    except KeyboardInterrupt:
	observer.stop()

    observer.join()