import functools
import logging

from io import BytesIO

from flask import Flask, jsonify, request, send_file

from config import settings
from core.images import render_image

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(name)s:%(lineno)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

app = Flask(__name__)


def response(data, status_code=200):
    return jsonify(data), status_code


def response_error(err_msg, status_code):
    return response({"error": err_msg}, status_code)


def required_arguments(*expected_args):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            actual_args = request.args
            for expected_arg in expected_args:
                if expected_arg not in actual_args:
                    return response_error(f"{expected_arg} name must be provided", 400)
                if not actual_args[expected_arg]:
                    return response_error(f"{expected_arg} is empty", 400)
            return func(*args, **kwargs)
        return wrapper
    return decorator


def handle_exception(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        # flake8: noqa
        except Exception as e:  # pylint: disable=broad-except
            return response_error(str(e), 500)
    return wrapper


@app.route("/render_image")
@handle_exception
def get_render_image():
    """Generate image using indigo or rdkit.

        query: required. Inchi or smiles of the compound.
        key: optional. Smiles to highlight compound in the requested. Possible if query is inchi.
    """
    query = request.args.get("query")
    key = request.args.get("key")
    if not query:
        return response_error("Query should be provided", 400)

    image = render_image(query, key)
    return send_file(BytesIO(image), mimetype="image/svg+xml")


@app.route("/substructure", methods=["GET"])
@handle_exception
@required_arguments("smiles")
def substructure_search():
    smiles = request.args.get("smiles")
    #TODO add finder
    matches = [smiles]
    return response(matches)


if __name__ == "__main__":
    app.run(debug=settings.DEBUG_MODE, host=settings.HOST, port=settings.PORT)
